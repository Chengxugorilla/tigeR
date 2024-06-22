#' @title Compare the roc curves of different models.
#' @description plot all the roc curves on a single plot.
#' @param ROC a list of roc objects that you want to campare.
#' @param lg.label a character vector specifying the legend labels, which must be of the same length as the ROC list.
#' @param line.col a character vector sepcifying the line colors.
#' @export

compare_roc <- function(ROC,lg.label = c("NB","SVM","RF","CC","ADB","LGB","LGT"),
                        line.col = c("black", "#F89B9B", "#4DB867", "#C64D6A", "#FDCEBC", "green", "#ADD8E6")){
  d <-
    lapply(seq_along(ROC), function(i){
      dplyr::arrange(
        data.frame(
          TPR = ROC[[i]]$sensitivities,
          FPR = ROC[[i]]$specificities,
          color = paste0("Curve ",i)
        ),.data$TPR)
    })

  lbs <- unlist(
    lapply(seq_along(ROC), function(i){
      paste0(lg.label[i]," ",sprintf("%.3f", ROC[[i]]$auc))
    }))

  plt <- ggplot()
  for (i in seq_along(d)){
    plt <- plt +
      geom_line(data = d[[i]], aes(x = .data$FPR, y = .data$TPR,
                                   color = .data$color), linewidth = 1)
  }

  plt +
    xlim(1, 0) +
    labs(x = "specificity", y = "sensitivity") +
    scale_color_manual(values = line.col[seq_along(lbs)],
                       labels = lbs) +
    theme(plot.background = element_rect(fill = "transparent", color = "transparent"),
          panel.background = element_rect(fill = "transparent", color = "black", linewidth = 1.5),
          panel.grid = element_blank(),
          legend.position = c(0.7,0.18),
          legend.background = element_rect(fill = "transparent",color = "transparent"),
          legend.box.background = element_rect(fill = "transparent",color = "transparent"),
          legend.key = element_blank(),
          legend.key.height = unit(0,"mm"),
          legend.key.spacing.y = unit(0,"mm"),
          legend.title = element_blank(),
          legend.text = element_text(size=7),
          legend.byrow = TRUE,
          axis.title = element_text(face = "bold", size = 12, color = "black"),
          axis.text = element_text(face = "bold", size = 9, color = "black"),
          aspect.ratio = 1) +
    guides(color=guide_legend(ncol = 1))
}
