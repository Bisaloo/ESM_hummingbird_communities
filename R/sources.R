ggplot2::theme_set(ggplot2::theme_minimal())

patchlong = c(
  "CR"="Crown",
  "BA"="Back",
  "RU"="Rump",
  "R1"="Tail",
  "TH"="Throat",
  "BR"="Breast",
  "BE"="Belly",
  "RE"="Wing"
)

cleanup_phylo = function(multiphylo) {

  multiphylo = mclapply(multiphylo, function(tree) {
    tree$tip.label = gsub('_', ' ', tree$tip.label)
    return(tree)
  })

  multiphylo = mclapply(multiphylo, function(tree) {
    tree$tip.label = gsub("Hylocharis", "Amazilia", tree$tip.label, fixed = TRUE)
    return(tree)
  })

  return(multiphylo)
}

create_spaco_commu = function(species) {
  # The spacodiR functions (PIst and TAUst) has a quite strict matrix input
  # format so we need to follow it.

  spaco_commu = dcast(species, Spname ~ idComm, value.var = 'CountOfSpname')
  rownames(spaco_commu) = spaco_commu$Spname

  # Remove Metallura odomae because it was absent from Paris and Lyon
  # collections and it could not be measured.
  spaco_commu = spaco_commu[spaco_commu$Spname!='Metallura odomae',
                            colnames(spaco_commu)!='Spname']
  colnames(spaco_commu) = paste0('id', unique(species$idComm))
  spaco_commu[is.na(spaco_commu)] = 0

  return(spaco_commu)

}

compute_null_spaco = function(sp.plot, phy = NULL, sp.traits = NULL,
                              method = '1s', nrand = 1000, set_seed = TRUE) {
  # SpacodiR doesn't include a statistics test to evaluate significance of the
  # computed indices. We randomize the community matrix and assume this
  # represent the index value under the null hypothesis.

  if (set_seed)
    set.seed(1993)

  replicate(nrand, {
    shuffled = eval(parse(text = paste0('resamp.', method, '(sp.plot)')))
    # Different shuffling methods are presented in Hardy 2008
    # https://doi.org/10.1111/j.1365-2745.2008.01421.x
    if (is.null(phy)) {
      null_spaco = spacodi.calc(sp.plot = shuffled, sp.traits = sp.traits)
    }
    else {
      null_spaco = spacodi.calc(sp.plot = shuffled, phy = phy)
    }
    # return PIst /  TAUst (depending on the case)
    return(null_spaco[[3]])
  })
}

summarize_spaco = function(spaco_res) {
  # Return a unique Xst and p value from a distribution of Xst and p

  mean_Xst = mean(spaco_res[,1])

  p_neg = compute_p_multiple(spaco_res[,2])
  p_pos = compute_p_multiple(spaco_res[,3])

  return(c(mean_Xst, p_neg, p_pos))
}

compute_p_multiple = function(p_distribution) {

  # Formula from Lou Jost
  # http://loujost.com/Statistics%20and%20Physics/Significance%20Levels/CombiningPValues.htm

  p_distribution = unlist(p_distribution)

  k = prod(p_distribution)

  if (k == 0) {
    # The formula won't work if k = 0. But the answer is obvious in this case.
    return(0)
  }
  else {
    vec_sum = seq(0, length(p_distribution)-1)
    # log() is ln
    vec_sum = ((- log(k)) ^ vec_sum) / factorial(vec_sum)

    p = k * sum(vec_sum)

    return(p)
  }
}

compute_PIst = function(sp.plot, multiphy) {

  # Compute a distribution of PIst and p
  res = mclapply(multiphy, function(tree) {

    pist = spacodi.calc(sp.plot = sp.plot, phy = tree)$PIst

    nrand = 1000

    pist_null = compute_null_spaco(sp.plot = sp.plot, phy = tree,
                                   nrand = nrand, set_seed = T)

    # Compute p values for both pist>0 AND pist<0
    pvalue_neg = sum(pist_null<pist) / nrand
    pvalue_pos = sum(pist_null>pist) / nrand

    return(c(pist, pvalue_neg, pvalue_pos))
  })

  res_df = do.call(rbind.data.frame, res)

  # Summarize the distribution of PIst and p in a single value
  PIst_with_p = summarize_spaco(res_df)

  return(PIst_with_p)
}

compute_TAUst_with_p = function(sp.plot, sp.traits = NULL) {

  taust = spacodi.calc(sp.plot = sp.plot,
                       sp.traits = sp.traits)$TAUst

  nrand = 1000

  taust_null = compute_null_spaco(sp.plot = sp.plot,
                                  sp.traits = sp.traits,
                                  nrand = nrand, set_seed = T)

  p_neg = sum(taust_null<taust) / nrand
  p_pos = sum(taust_null>taust) / nrand

  return(c(taust, p_neg, p_pos))

}

compute_TAUst_onlist = function(sp.plot, liste_sp.trait) {

  taust_p = mclapply(liste_sp.trait, function(sp.trait) {
    compute_TAUst_with_p(sp.plot = sp.plot, sp.trait = sp.trait)
  })

  taust_p = do.call(rbind.data.frame, taust_p)
  colnames(taust_p) = c("TAUst", "p_neg", "p_pos")
  rownames(taust_p) = names(liste_sp.trait)

  return(taust_p)
}

compute_decoupled = function(list_phylo, todecouple) {

  # Create a new trait dataframe after "removing" the effect of the phylogeny
  # on the phenotypic structure.

  # Be careful, todecouple HAS TO be a dataframe, a named matrix won't cut it.
  if (!is.data.frame(todecouple))
    stop("todecouple must be a dataframe")

  decoupled = mclapply(list_phylo, function(tree) {

    pruned_tree = drop.tip(tree, tree$tip.label[!tree$tip.label %in% rownames(todecouple)])

    decouple(traits = todecouple[match(pruned_tree$tip.label,
                                 rownames(todecouple)),,drop=FALSE],
             tree = pruned_tree,
             keep_axes_method = "lm")$dcFdist
  })

  # Because we decouple the traits from the phylogeny by using a tree set
  # (instead of a consensus), we end up with a list (size = number of trees)
  # where each element is a dataframe of decoupled traits.
  return(decoupled)

}

compute_decoupled_onlist = function(list_phylo, list_todecouple){

  lapply(list_todecouple, function(df_todecouple) {
    compute_decoupled(list_phylo = list_phylo, todecouple = df_todecouple)
  })

}

spaco_decoupled = function(sp.plot, sp.decoupled) {

  # Loop over the list of decoupled df (size = number of trees)

  res_eachtree = mclapply(sp.decoupled, function(df_decoupled) {
    compute_TAUst_with_p(sp.plot = sp.plot, sp.traits = df_decoupled)})

  mat_eachtree = do.call(rbind.data.frame, res_eachtree)

  summarize_spaco(mat_eachtree)

}

spaco_decoupled_onlist = function(sp.plot, liste_sp.trait) {

  # Same as previously but with nested list: by area > by tree > decoupled df

  dctaust_p = lapply(liste_sp.trait,
                   function(sp.trait) {
    spaco_decoupled(sp.plot = sp.plot, sp.decoupled = sp.trait)
  })

  dctaust_p = do.call(rbind.data.frame, dctaust_p)
  colnames(dctaust_p) = c("dcTAUst", "p_neg", "p_pos")
  rownames(dctaust_p) = names(liste_sp.trait)

  return(dctaust_p)
}

compute_spectra_all = function(data_all) {

  # Parse data into a tidy rspec object

  spectra_df = data.frame(seq(300,700),
                          t(data_all[,paste0("X", seq(300,700))]))
  colnames(spectra_df) = c("wl",
                           paste(data_all$Species, data_all$Patch, data_all$Exp, sep = "_"))
  rownames(spectra_df) = NULL

  spectra_df = pavo::as.rspec(spectra_df)
  spectra_df = pavo::procspec(spectra_df, opt = "smooth", fixneg = "zero")

  # Only keep one subspecies when more than one was measured
  spectra_df = spectra_df[, !grepl("hesperus|torquata_torquata|jamesonii|leucura|albogularis",
                                   colnames(spectra_df))]

  return(spectra_df)
}

compute_spectra_max = function(spectra_df) {

  spectra_max = spectra_df[, grepl("(wl|MAX|NI)$", colnames(spectra_df))]

  colnames(spectra_max) = sapply(strsplit(colnames(spectra_max), "_"),
                                 function(x) paste(x[c(1,2,length(x)-1)], collapse = "_"))

  colnames(spectra_max)[1] = "wl"

  return(spectra_max)

}

compute_brightness_max = function(tetra_all) {

  brightness_max = tetra_all[tetra_all$Exp %in% c("NI", "MAX"),
                             c("Species", "Patch", "lum")]


  return(brightness_max)
}

compute_vismodel_all = function(spectra_all, habitat_df, vision_df, illum_C, illum_U) {

  visdf = vision_df[,c("wl","u","s","m","l")]
  achdf = vision_df[,"dc"]

  liste_species_spectra_all = sapply(strsplit(colnames(spectra_all), "_"),
                                     function(x) paste(x[c(1,2)], collapse = " "))

  spectra_C = spectra_all[,liste_species_spectra_all %in% habitat_df$TrueName[habitat_df$SimpleStrata=="C"]]
  spectra_C$wl = spectra_all$wl

  spectra_U = spectra_all[,liste_species_spectra_all %in% habitat_df$TrueName[habitat_df$SimpleStrata=="U"]]
  spectra_U$wl = spectra_all$wl

  vm_C = vismodel(spectra_C, qcatch = "Qi", visual = visdf, achromatic = achdf,
                  illum = illum_C$largegap, relative = TRUE)

  vm_U = vismodel(spectra_U, qcatch = "Qi", visual = visdf, achromatic = achdf,
                  illum = illum_U$LuxForestShade, relative = TRUE)

  vm_all = rbind(vm_C, vm_U)
}

compute_vismodel_max = function(vm_all) {
  vm_max = vm_all[grepl("_MAX$", rownames(vm_all)),]

  vm_max$Species = gsub("^([[:alpha:]]+_[[:alpha:]]+).*", "\\1", rownames(vm_max))
  vm_max$Patch   = gsub(".*_([[:upper:]|[:digit:]]+)_MAX$", "\\1", rownames(vm_max))

  return(vm_max)
}

compute_tetra_all = function(vm_all) {

  tetra_all = colspace(vismodeldata = vm_all, space = "tcs")

  tetra_all$Species = sapply(strsplit(rownames(tetra_all), "_"),
                             function(x) paste(x[c(1,2)], collapse = " "))
  tetra_all$Patch = sapply(strsplit(rownames(tetra_all), "_"),
                           function(x) x[length(x)-1])
  tetra_all$Exp = sapply(strsplit(rownames(tetra_all), "_"),
                         function(x) x[length(x)])

  rownames(tetra_all) = NULL

  return(tetra_all)
}

compute_contrasts = function(spectra, habitat_df, vision_df, illum_C, illum_U, bg_C, bg_U) {

  visdf = vision_df[,c("wl","u","s","m","l")]
  achdf = vision_df[,"dc"]

  liste_species_spectra = sapply(strsplit(colnames(spectra), "_"),
                                 function(x) paste(x[c(1,2)], collapse = " "))

  spectra_C = spectra[,liste_species_spectra %in% habitat_df$TrueName[habitat_df$SimpleStrata=="C"]]
  spectra_C$wl = spectra$wl

  spectra_U = spectra[,liste_species_spectra %in% habitat_df$TrueName[habitat_df$SimpleStrata=="U"]]
  spectra_U$wl = spectra$wl

  bg_C = procspec(as.rspec(bg_C),
                  opt = I("smooth"),
                  fixneg = I("zero"))
  bg_U = procspec(as.rspec(bg_U),
                  opt = I("smooth"),
                  fixneg = I("zero"))

  df_C = cbind(spectra_C, bg_C)
  df_U = cbind(spectra_U, bg_U)

  vm_C = vismodel(df_C, qcatch = "fi", visual = visdf, achromatic = achdf,
                  illum = illum_C$largegap, relative = TRUE)

  vm_U = vismodel(df_U, qcatch = "fi", visual = visdf, achromatic = achdf,
                  illum = illum_U$LuxForestShade, relative = TRUE)

  vm_all = rbind(vm_C, vm_U)

  # From Hart 2000
  contrasts = coldist(vm_all, noise = 'neural', achro = TRUE, n = c(1, 1.9, 2.2, 2.1))

  return(contrasts)
}

prepare_tetra_max = function(tetra_all) {

  # We consider that we get the "true" colour for the angle configuration where
  # the brightness is maximal.
  # Because we will use this a lot, it is better to keep it in a separate
  # dataframe.

  tetra_max = tetra_all[tetra_all$Exp %in% c("NI", "MAX"),
                        c("Species", "Patch", "x", "y", "z")]

  return(tetra_max)
}

compute_crossvolmat = function(tetra_max) {

  # Returns a distance matrix taking species as rows and columns with each cell
  # indicating the volume overlap (in the tetrahedron) of the two species colour
  # volumes.

  species_list = unique(tetra_max$Species)

  all_crossvol = combn(species_list, m = 2, FUN = function(pair) {
    tryCatch(
    voloverlap(tetra_max[tetra_max$Species == pair[1], c("x", "y", "z")],
               tetra_max[tetra_max$Species == pair[2], c("x", "y", "z")])$overlapvol,
    error = function(e) NA)
    })

  crossvolmat = matrix(NA, nrow = length(species_list),
                           ncol = length(species_list),
                       dimnames = rep(list(species_list),2))

  crossvolmat[upper.tri(crossvolmat)] = all_crossvol
  crossvolmat[lower.tri(crossvolmat)] = all_crossvol

  diag(crossvolmat) = by(tetra_max, tetra_max$Species, function(tetra_species) {
    ifelse(nrow(tetra_species)>3,
      convhulln(tetra_species[, c("x", "y", "z")], "FA")$vol,
      NA
    )
  })

  return(crossvolmat)
}

plot_comcolvol = function(tetra_max, spaco_commu) {

  total_colvol = convhulln(tetra_max[, c("x", "y", "z")], "FA")$vol

  res = apply(X = spaco_commu,
        MARGIN = 2,
        FUN = function(commu_compo) {
          spcom_i = rownames(spaco_commu)[as.logical(commu_compo)]
          tetracom_i = tetra_max[tetra_max$Species %in% spcom_i,]
          colvol_i = convhulln(tetracom_i[, c("x", "y", "z")], "FA")$vol
          nbsp_i = sum(commu_compo)

          return(c(nbsp_i, colvol_i))
        })

  res = data.frame(t(res))
  colnames(res) = c("nbsp", "colvol")
#
#   m_comcolvol_lin = lm(colvol ~ nbsp, res)
#   m_comcolvol_log = lm(colvol ~ log(nbsp), res)
#
#   AIC_lin = AIC(m_comcolvol_lin)
#   AIC_log = AIC(m_comcolvol_log)
#
#   r_lin = summary(m_comcolvol_lin)$adj.r.squared
#   r_log = summary(m_comcolvol_log)$adj.r.squared
#
#   res$lin = m_comcolvol_lin$fitted.values
#   res$log = m_comcolvol_log$fitted.values
#
#   ggplot(res, aes(x = nbsp, y = colvol)) +
#     geom_point() +
#     xlab('Species count in community') + ylab('colour volume of community') +
#     geom_line(aes(y = lin, linetype = "Linear")) +
#     annotate("text", label = sprintf("R² = %.2f, AIC = %.2f", r_lin, AIC_lin), x = max(res$nbsp)-4, y = max(res$lin)+0.001) +
#     geom_line(aes(y = log, linetype = "Log")) +
#     annotate("text", label = sprintf("R² = %.2f, AIC = %.2f", r_log, AIC_log), x = max(res$nbsp)-4, y = max(res$log)-0.006) +
#     scale_colour_brewer("", palette="Set1") +
#     geom_hline(yintercept = total_colvol) +
#     annotate("text", label='colour volume of all communities pooled',  x = median(res$nbsp), y = total_colvol+0.002) +
#     xlim(c(0, max(res$nbsp))) +
#     ylim(c(0, total_colvol+0.005)) +
#     theme(legend.position = "none")

  nullcolvol = lapply(seq(1, max(res$nbsp)), function(nbsp) {
    replicate(10000, {
    null_sp = sample(unique(tetra_max$Species), size = nbsp)
    null_tcs = tetra_max[tetra_max$Species %in% null_sp, ]
    null_vol = convhulln(null_tcs[, c("x", "y", "z")], "FA")$vol
    })
  })

  nullcolvol = do.call(rbind.data.frame, nullcolvol)
  rownames(nullcolvol) = seq(1, max(res$nbsp))

  nullres = data.frame(t(apply(nullcolvol, 1, function(x)
    c("colvol" = mean(x), "avghi" = mean(x) + 2*sd(x),
                       "avglo" = mean(x) - 2*sd(x)))))
  nullres$nbsp = seq(1, max(res$nbsp))

  ggplot(res, aes(x = nbsp, y = colvol)) +
  geom_ribbon(data = nullres,
              aes(ymin = avghi, ymax = avglo), fill = "grey80", alpha = 0.5) +
  geom_line(data = nullres) +
  geom_point() +
  xlab('Species count in community') +
  ylab('Colour volume of community') +
  labs(tag = "(a)")

}

plot_comm_brightnessrange = function(brightness_max, spaco_commu) {

  res = apply(X = spaco_commu,
              MARGIN = 2,
              FUN = function(commu_compo) {
                spcom_i = rownames(spaco_commu)[as.logical(commu_compo)]
                brightness_i = brightness_max$lum[brightness_max$Species %in% spcom_i]
                max_brightness_i = max(brightness_i)
                min_brightness_i = min(brightness_i)
                range_brightness_i = max_brightness_i - min_brightness_i
                nbsp_i = sum(commu_compo)

                return(c(nbsp_i, max_brightness_i, range_brightness_i))
              })

  res = data.frame(t(res))
  colnames(res) = c("nbsp", "max_brightness", "range_brightness")

  # m_comm_bright_lin = lm(range_brightness ~ nbsp, res)
  # m_comm_bright_log = lm(range_brightness ~ log(nbsp), res)
  #
  # AIC_lin = AIC(m_comm_bright_lin)
  # AIC_log = AIC(m_comm_bright_log)
  #
  # r_lin = summary(m_comm_bright_lin)$adj.r.squared
  # r_log = summary(m_comm_bright_log)$adj.r.squared
  #
  # res$lin = m_comm_bright_lin$fitted.values
  # res$log = m_comm_bright_log$fitted.values
  #
  # ggplot(res, aes(x = nbsp, y = range_brightness)) +
  #   geom_point() +
  #   xlab("Species count in community") +
  #   ylab("Range of brightness in the community") +
  #   geom_line(aes(y = lin, linetype = "Linear")) +
  #   annotate("text", label = sprintf("R² = %.2f, AIC = %.2f", r_lin, AIC_lin), x = max(res$nbsp)-4, y = max(res$lin)+50) +
  #   geom_line(aes(y = log, linetype = "Log")) +
  #   annotate("text", label = sprintf("R² = %.2f, AIC = %.2f", r_log, AIC_log), x = max(res$nbsp)-4, y = max(res$log)-160) +
  #   scale_colour_brewer("", palette="Set1") +
  #   xlim(c(0, max(res$nbsp))) +
  #   theme(legend.position = "none")

  nullcomb = lapply(seq(1, max(res$nbsp)), function(nbsp) {
    replicate(10000, {
      null_sp = sample(unique(brightness_max$Species), size = nbsp)
      null_b = brightness_max$lum[brightness_max$Species %in% null_sp]
      max_null_b = max(null_b)
      min_null_b = min(null_b)
      range_null_b = max_null_b - min_null_b
      return(range_null_b)
    })
  })

  nullcomb = do.call(rbind.data.frame, nullcomb)
  rownames(nullcomb) = seq(1, max(res$nbsp))

  nullres = data.frame(t(apply(nullcomb, 1, function(x)
    c("range_brightness" = mean(x), "avghi" = mean(x) + 2*sd(x),
      "avglo" = mean(x) - 2*sd(x)))))
  nullres$nbsp = seq(1, max(res$nbsp))

  ggplot(res, aes(x = nbsp, y = range_brightness)) +
    geom_ribbon(data = nullres,
                aes(ymin = avghi, ymax = avglo), fill = "grey80", alpha = 0.5) +
    geom_line(data = nullres) +
    geom_point() +
    xlab('Species count in community') +
    ylab('Brightness range of community') +
    labs(tag = "(b)")

}

create_area_patch_sp = function(area_patches, corpatsp) {

  area_grouped_patches = apply(corpatsp, 2, function(x) {
    pattern_bylines = split(names(x), x)

    group_patches = lapply(pattern_bylines, unique)

    area_group_patches = sapply(group_patches, function(group) {
      area_group = sum(area_patches$Area[area_patches$Patch %in% group])
    })

    return(area_group_patches)
  })

  return(area_grouped_patches)
}

create_bary = function(tetra_max, area_patches, corpatsp) {

  # Compute the barycenter of the colour volume for each species by weighting
  # by the area of the patch where this colour was measured.

  # TODO: this only take into account the measured patches. It should also use
  # the patterns to simulate non-measured colours.

  area_grouped_patches = create_area_patch_sp(area_patches, corpatsp)

  list_sp = unique(tetra_max$Species)

  bary_all = lapply(list_sp, function(species) {
    tetra_species = tetra_max[tetra_max$Species == species,]
    area_patches_species = get(species, area_grouped_patches)
    area_patches_species = data.frame("Patch" = names(area_patches_species),
                                      "Area" = area_patches_species,
                                      row.names = NULL)

    tetra_area_species = merge(tetra_species, area_patches_species, by = "Patch")

    bary = colSums(tetra_area_species$Area * tetra_area_species[,c("x","y","z")]) / sum(tetra_area_species$Area)

    return(bary)
  })

  bary_all = data.frame(do.call(rbind, bary_all))

  rownames(bary_all) = list_sp

  return(bary_all)
}

compute_dist_bary = function(barycenters) {
  dist(x = barycenters, method = 'euclidean', diag = TRUE, upper = TRUE)
}

partial_mantel = function(list_phylo, distmat_1, distmat_2) {

  distmat_1 = as.matrix(distmat_1)
  distmat_2 = as.matrix(distmat_2)

  sp_list = sort(intersect(rownames(distmat_1), rownames(distmat_2)))

  # Make sure the species are in the same order in both matrices
  distmat_1 = distmat_1[match(sp_list, rownames(distmat_1)),
                        match(sp_list, colnames(distmat_1))]
  distmat_2 = distmat_2[match(sp_list, rownames(distmat_2)),
                        match(sp_list, colnames(distmat_2))]
  list_phylo = mclapply(list_phylo, function(tree) {
    drop.tip(tree, tree$tip.label[!tree$tip.label %in% sp_list])
  })

  # Do partial mantel test on a tree distribution and return a single summary
  # p-value
  # vegan::mantel.partial already use parallel processing
  res_pmantel = lapply(list_phylo, function(tree){
    mantel.partial(cophenetic(tree),
                   distmat_1,
                   distmat_2)
  })
  R_pmantel = mean(sapply(res_pmantel, function(x) x$statistic))
  p_pmantel = compute_p_multiple(sapply(res_pmantel, function(x) x$signif))

  res_mantel = vegan::mantel(distmat_1, distmat_2)
  R_mantel = res_mantel$statistic
  p_mantel = res_mantel$signif

  return(c("R_mantel" = R_mantel,
           "p_mantel" = p_mantel,
           "R_pmantel" = R_pmantel,
           "p_pmantel" = p_pmantel))
}

partial_mantel_onlist = function(list_phylo, distmat_1, list_distmat_2) {

  pmantel_list = lapply(list_distmat_2, function(distmat_i) {
    partial_mantel(list_phylo, distmat_1, distmat_i)
  })

  pmantel_res = do.call(rbind.data.frame, pmantel_list)
  colnames(pmantel_res) = c("R_mantel", "p_mantel", "R_pmantel", "p_pmantel")
  rownames(pmantel_res) = names(list_distmat_2)

  return(pmantel_res)
}

create_Cdist = function(spaco_commu) {
  Cdist = vegan::designdist(spaco_commu, method = "(A-J)*(B-J)/(A+B-J)", terms = "binary")

  return(Cdist)
}

colvol_commu_area = function(tetra_area, spaco_commu) {

  res = t(apply(X = spaco_commu,
              MARGIN = 2,
              FUN = function(commu_compo) {
                spcom_i = rownames(spaco_commu)[as.logical(commu_compo)]
                tetracom_i = tetra_area[rownames(tetra_area) %in% spcom_i,]
                colvol_i = ifelse(nrow(tetracom_i)>3,
                                  yes = convhulln(tetracom_i[, c("x","y","z")], "FA")$vol,
                                  no = NA)
                nbsp_i = sum(commu_compo)

                return(c(nbsp_i, colvol_i))
              }))

  colnames(res) = c('nbsp', 'colvol')

  return(res)
}

plot_commu_colvol_patch = function(tetra_byarea, spaco_commu) {

  colvol_commu_patches = mclapply(tetra_byarea, function(tetra_area) {
    colvol_commu_area(tetra_area, spaco_commu)})

  colvol_commu_patches = do.call(rbind.data.frame, colvol_commu_patches)
  colvol_commu_patches$patch = sapply(strsplit(rownames(colvol_commu_patches), "\\."),
                                      function(x) x[1])

  ggplot(colvol_commu_patches, aes(x = nbsp, y = colvol)) + geom_point() +
    facet_wrap( ~ patch) +
    scale_colour_brewer(palette = "Set1") +
    xlab("Species count in community") +
    ylab("Mean specific colour volume")
}

compute_huedist_patch = function(tetra_byarea) {

  dist_patch = mclapply(tetra_byarea, function(tetra_area) {
    dist(tetra_area, method = "euclidean")})

  return(dist_patch)
}

prepare_IR_MAX = function(tetra_all) {

  IR_MAX = do.call(rbind.data.frame, by(tetra_all, paste(tetra_all$Species, tetra_all$Patch, sep='_'), function(tetra_sppatch) {
    MAX     = tetra_sppatch[tetra_sppatch$Exp %in% c("MAX", "NI"), c("x", "y", "z")]
    IR      = tetra_sppatch[tetra_sppatch$Exp %in% c("IR2", "NI"), c("x", "y", "z")]
    Patch   = unique(tetra_sppatch[, "Patch"])
    Species = unique(tetra_sppatch[, "Species"])

    return(c(MAX, IR, Patch, Species))
  }))

  colnames(IR_MAX) = c('MAX_x', 'MAX_y', 'MAX_z',
                       'IR_x', 'IR_y', 'IR_z',
                       "Patch",
                       "Species")

  return(IR_MAX)
}

prepare_iri = function(tetra_all) {

  iri = do.call(rbind.data.frame, by(tetra_all, paste(tetra_all$Species, tetra_all$Patch, sep='_'), function(tetra_sppatch) {
    MAX     = tetra_sppatch[tetra_sppatch$Exp %in% c("MAX", "NI"), c("x", "y", "z")]
    IR      = tetra_sppatch[tetra_sppatch$Exp %in% c("IR2", "NI"), c("x", "y", "z")]
    IRI     = MAX - IR
    Patch   = unique(tetra_sppatch[, "Patch"])
    Species = unique(tetra_sppatch[, "Species"])

    return(c(MAX, IRI, Patch, Species))
  }))

  colnames(iri) = c('MAX_x', 'MAX_y', 'MAX_z',
                       'IRI_x', 'IRI_y', 'IRI_z',
                       "Patch",
                       "Species")

  return(iri)
}

prepare_df_byarea = function(list_areas, df) {

  df_byarea = lapply(list_areas, function(area) {
    df_area = df[df$Patch==area,]
    rownames(df_area) = df_area$Species

    df_area = df_area[, !colnames(df_area) %in% c("Species", "Patch"), drop = FALSE]

    return(df_area)
  })

  names(df_byarea) = list_areas

  return(df_byarea)
}

compute_cooc = function(spaco_commu) {

  res = data.frame("Species" = rownames(spaco_commu),
                   "nbcommu" = rowSums(spaco_commu),  # number of communities where this species has been reported
                   "nbcooc" = rowSums((as.matrix(spaco_commu) %*% t(as.matrix(spaco_commu)))!=0) # number of different species co-ocurring with the focal species
                  )
  res$nbcoocmoy = res$nbcooc / res$nbcommu

  return(res)
}

plot_cooc_colvol = function(crossvolmat, species_cooc) {

  species_cooc$colvol = diag(crossvolmat)

  ggplot(data = species_cooc, aes(x = nbcoocmoy, y = colvol)) +
    geom_point()
}

plot_map = function(vector, raster, communities) {

  raster_spdf = as(raster, "SpatialPixelsDataFrame")
  raster_df = as.data.frame(raster_spdf)
  colnames(raster_df) = c("value", "x", "y")

  ggplot() +
    geom_tile_rast(data=raster_df, aes(x=x, y=y, fill=value)) +
    scale_fill_viridis(name = "Elevation") +
    geom_polygon(data = vector, aes(long, lat, group=group), fill = 'transparent') +
    geom_point(data = communities, aes(LongDecDeg, LatDecDeg), color = 'red') +
    coord_equal() + xlab('Longitude') + ylab('Latitude') +
    theme(legend.position = 'bottom') +
    xlim(range(communities$LongDecDeg) + c(-0.5, 0.5)) +
    ylim(range(communities$LatDecDeg) + c(-0.5, 0.5)) +
    annotate('text', x = -76.1, y = -5, label = '100km') +
    annotate('segment', x = -76.65, xend = -75.55, y = - 4.75, yend = -4.75)
}

plot_coverage = function(tree, spaco_commu) {

  if (class(tree)=="list") {
    tree = ape::consensus(tree)
  }
  tree = drop.tip(tree, c('Hemiprocne_mystacea', 'Aegotheles_insignis'))
  tree = chronoMPL(tree)

  list_species = rownames(spaco_commu)

  tree$tip.label = gsub("[[:space:]]", "_", tree$tip.label)
  list_species = gsub("[[:space:]]", "_", list_species)

  group_struct = list(
    'ABSENT' = c(),
    'PRESENT' = list_species
  )
  tree = groupOTU(tree, group_struct)

  ggtree(tree, layout="circular", aes(color = group)) +
    geom_tiplab2(size = 1) +
    scale_color_manual(values = c("grey", "firebrick")) +
    geom_treescale(x = 42, fontsize = 2)
}

import_data_all = function() {

  data_all = read.csv(file = "data/Data_all.csv")
#  subset(data_all, Angle.separating.fibers <= 11)
}


plot_col_tree = function(spectra_max, tree) {

  consensus = maxCladeCred(tree)

  p = ggtree(consensus, ladderize = FALSE) + theme_tree2()

  dd <- spectra_max %>%
    summary(subset = "H1") %>%
    rownames_to_column() %>%
    filter(grepl("_(CR|BA|RU|R1|TH|BR|BE|RE)$", rowname)) %>%
    separate(rowname, "_", into = c("genus", "species", "patch")) %>%
    mutate(taxa = paste(genus, species)) %>%
    rowid_to_column() %>%
    dplyr::select(taxa, patch, rowid) %>%
    mutate(taxa = factor(taxa, levels = consensus$tip.label),
           patch = patchlong[patch])

  colp = ggplot(dd, aes(x = patch, y = taxa, color= factor(rowid))) + geom_point(shape = "square", size = 2) + theme(legend.position = "none") +
    scale_color_manual(values = unname(spec2rgb(spectra_max[, c(1, grep("_(CR|BA|RU|R1|TH|BR|BE|RE)$", colnames(spectra_max)))]))) +
    ylab("") +
    theme(axis.text.y = element_text(hjust = 0.5))

  multiplot(p, colp, ncol = 2)
}
