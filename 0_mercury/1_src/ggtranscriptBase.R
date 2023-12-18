
plot_simply <- function(merged){
  gexons = merged
  gintrons = gexons %>% to_intron(group_var = "transcript_id")
  grescaled = shorten_gaps(gexons, gintrons, group_var = "transcript_id")
  
  p <- grescaled %>%
    dplyr::filter(type == "exon") %>%
    ggplot(aes(xstart = start,xend = end,y = transcript_id)) +
    geom_range(aes(fill = Type)) +
    labs(y ="") + 
    geom_intron(data = grescaled %>% dplyr::filter(type == "intron"),arrow.min.intron.length = 300) +
    theme_classic() +
    theme(legend.position = "None", 
          axis.line.x = element_line(colour = "grey80"),
          panel.background = element_rect(fill = "white", colour = "grey50"),
          panel.border = element_rect(fill = NA, color = "grey50", linetype = "dotted"),
          axis.text.y= element_text(size=12),
          strip.text.y = element_text(size = 12, color = "black"),
          strip.background = element_rect(fill = "white", colour = "grey50")) + scale_fill_manual(values = c("red","black"))
  
  return(p)
}

geom_range <- function(mapping = NULL, data = NULL,
                       stat = "identity", position = "identity",
                       ...,
                       vjust = NULL,
                       linejoin = "mitre",
                       na.rm = FALSE,
                       show.legend = NA,
                       inherit.aes = TRUE) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomRange,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      vjust = vjust,
      linejoin = linejoin,
      na.rm = na.rm,
      ...
    )
  )
}

geom_intron <- function(mapping = NULL, data = NULL,
                        stat = "identity", position = "identity",
                        ...,
                        arrow = grid::arrow(ends = "last", length = grid::unit(0.1, "inches")),
                        arrow.fill = NULL,
                        lineend = "butt",
                        linejoin = "round",
                        na.rm = FALSE,
                        arrow.min.intron.length = 0,
                        show.legend = NA,
                        inherit.aes = TRUE) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomIntron,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      arrow = arrow,
      arrow.fill = arrow.fill,
      lineend = lineend,
      linejoin = linejoin,
      na.rm = na.rm,
      arrow.min.intron.length = arrow.min.intron.length,
      ...
    )
  )
}

#' `GeomIntron` is pretty much `ggplot2::GeomSegment` with the `required_aes`
#' changed to `xstart`/`xend` to match genetic nomenclature and the added arrows
#' to indicate direction of transcription (configured with `strand` and
#' `arrow.min.intron.length`)
#' @noRd
GeomIntron <- ggplot2::ggproto("GeomIntron", ggplot2::GeomSegment,
                               required_aes = c("xstart", "xend", "y"),
                               default_aes = aes(
                                 colour = "black",
                                 size = 0.5,
                                 linetype = 1,
                                 alpha = NA,
                                 strand = "+"
                               ),
                               setup_params = function(data, params) {
                                 # check that arrow.min.intron.length numeric is >= 0
                                 arrow.min_numeric <- is.numeric(params$arrow.min.intron.length)
                                 arrow.min_neg <- params$arrow.min.intron.length < 0
                                 
                                 if (!arrow.min_numeric | arrow.min_neg) {
                                   stop("arrow.min.intron.length must be a numeric > 0")
                                 }
                                 
                                 params
                               },
                               setup_data = function(data, params) {
                                 # needed to permit usage of xstart/xend
                                 transform(
                                   data,
                                   x = xstart,
                                   yend = y,
                                   xstart = NULL
                                 )
                               },
                               draw_panel = function(data,
                                                     panel_params,
                                                     coord,
                                                     arrow = NULL,
                                                     arrow.fill = NULL,
                                                     lineend = "butt",
                                                     linejoin = "round",
                                                     na.rm = FALSE,
                                                     arrow.min.intron.length = 0) {
                                 
                                 # check that strand is scalar and one of "+" or "-"
                                 .check_strand(data$strand)
                                 
                                 # first, create the intron grob, which is just a pure line (no arrow)
                                 intron_grob <- ggplot2::GeomSegment$draw_panel(
                                   data = data,
                                   panel_params = panel_params,
                                   coord = coord,
                                   arrow = NULL,
                                   arrow.fill = NULL,
                                   lineend = lineend,
                                   linejoin = linejoin,
                                   na.rm = na.rm
                                 )
                                 
                                 # then, create the arrow grobs, one per strand
                                 # need both as the direction of arrow (as far I can tell) is
                                 # is dependent on the orientation of the x/xend
                                 strand_arrow_plus_grob <- .create_strand_arrow_grob(
                                   target_strand = "+",
                                   arrow.min.intron.length = arrow.min.intron.length,
                                   data = data,
                                   panel_params = panel_params,
                                   coord = coord,
                                   arrow = arrow,
                                   arrow.fill = arrow.fill,
                                   lineend = lineend,
                                   linejoin = linejoin,
                                   na.rm = na.rm
                                 )
                                 
                                 strand_arrow_minus_grob <- .create_strand_arrow_grob(
                                   target_strand = "-",
                                   arrow.min.intron.length = arrow.min.intron.length,
                                   data = data,
                                   panel_params = panel_params,
                                   coord = coord,
                                   arrow = arrow,
                                   arrow.fill = arrow.fill,
                                   lineend = lineend,
                                   linejoin = linejoin,
                                   na.rm = na.rm
                                 )
                                 
                                 # draw_panel expects return of a grob
                                 # here, as we build multiple grobs (i.e. intron lines + arrows)
                                 # we use a grobTree to combine the two
                                 grid::grobTree(
                                   intron_grob,
                                   strand_arrow_plus_grob,
                                   strand_arrow_minus_grob
                                 )
                               }
)

#' @keywords internal
#' @noRd
.check_strand <- function(strand) {
  # TODO - add option for "*" arrow?
  any_na <- any(is.na(strand))
  plus_minus <- !(all(strand %in% c("+", "-")))
  
  if (any_na | plus_minus) {
    stop("strand values must be one of '+' and '-'")
  }
  
  return(invisible())
}

to_intron <- function(exons, group_var = NULL) {
  .check_coord_object(exons)
  .check_group_var(exons, group_var)
  
  # TODO - switch this to using GenomicRanges::gaps()?
  
  if (!is.null(group_var)) {
    exons <- exons %>% dplyr::group_by_at(.vars = group_var)
  }
  
  # make sure exons are arranged by coord, so that dplyr::lag works correctly
  exons <- exons %>%
    dplyr::arrange(start, end)
  
  # obtain intron start and ends
  introns <- exons %>%
    dplyr::mutate(
      intron_start := dplyr::lag(end),
      intron_end := start,
      type = "intron"
    ) %>%
    dplyr::select(-start, -end)
  
  # remove the introduced artifact NAs
  introns <- introns %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(intron_start) & !is.na(intron_end))
  
  # filter out introns with a width of 1, this should only happen when
  # utrs are included and are directly adjacent to end of cds
  introns <- introns %>% dplyr::filter(abs(intron_end - intron_start) != 1)
  
  introns <- introns %>% dplyr::rename(start = intron_start, end = intron_end)
  
  return(introns)
}

#' @keywords internal
#' @noRd
.create_strand_arrow_grob <- function(target_strand,
                                      arrow.min.intron.length,
                                      data,
                                      panel_params,
                                      coord,
                                      arrow,
                                      arrow.fill,
                                      lineend,
                                      linejoin,
                                      na.rm) {
  
  # filter for introns that match target strand
  # and have a length above arrow.min.intron.length
  match_strand <- data$strand == target_strand
  ab_min <- abs(data$x - data$xend) > arrow.min.intron.length
  arrow_data <- data[match_strand & ab_min, ]
  
  # if there are no arrows to plot, use a nullGrob() to add nothing
  if (nrow(arrow_data) == 0) {
    arrow_grob <- grid::nullGrob()
  } else {
    
    # obtain the the correct orientation of arrow (dependent on strand)
    # as the arrow can only be placed at either end of a geom_segment/path
    # the strand changes the x/xends around, shifting the around direction
    if (target_strand == "+") {
      arrow_data <- transform(
        arrow_data,
        xend = (x + xend) / 2
      )
    } else {
      arrow_data <- transform(
        arrow_data,
        mid = (x + xend) / 2,
        x = xend
      )
      arrow_data <- transform(
        arrow_data,
        xend = mid
      )
    }
    
    arrow_grob <- ggplot2::GeomSegment$draw_panel(
      data = arrow_data,
      panel_params = panel_params,
      coord = coord,
      arrow = arrow,
      arrow.fill = arrow.fill,
      lineend = lineend,
      linejoin = linejoin,
      na.rm = na.rm
    )
  }
  
  return(arrow_grob)
}

GeomRange <- ggplot2::ggproto("GeomRange", ggplot2::GeomTile,
                              required_aes = c("xstart", "xend", "y"),
                              default_aes = aes(
                                fill = "grey",
                                colour = "black",
                                size = 0.25,
                                linetype = 1,
                                alpha = NA,
                                height = NA
                              ),
                              setup_data = function(data, params) {
                                # modified from ggplot2::GeomTile
                                #data$height <- data$height %||% params$height %||% 0.5
                                data$height <- ifelse(is.null(data$height), ifelse(is.null(params$height), 0.5, params$height), data$height)
                                
                                
                                transform(
                                  data,
                                  xmin = xstart,
                                  xmax = xend,
                                  ymin = y - height / 2,
                                  ymax = y + height / 2,
                                  height = NULL
                                )
                              },
                              draw_panel = function(self,
                                                    data,
                                                    panel_params,
                                                    coord,
                                                    vjust = NULL,
                                                    lineend = "butt",
                                                    linejoin = "mitre") {
                                if (!coord$is_linear()) {
                                  # prefer to match geom_curve and warn
                                  # rather than copy the implementation from GeomRect for simplicity
                                  # also don'think geom_range would be used for non-linear coords
                                  warn("geom_ is not implemented for non-linear coordinates")
                                }
                                
                                coords <- coord$transform(data, panel_params)
                                grid::rectGrob(
                                  coords$xmin, coords$ymax,
                                  width = coords$xmax - coords$xmin,
                                  height = coords$ymax - coords$ymin,
                                  default.units = "native",
                                  just = c("left", "top"),
                                  vjust = vjust,
                                  gp = grid::gpar(
                                    col = coords$colour,
                                    fill = ggplot2::alpha(coords$fill, coords$alpha),
                                    lwd = coords$size * ggplot2::.pt,
                                    lty = coords$linetype,
                                    linejoin = linejoin,
                                    lineend = lineend
                                  )
                                )
                              }
)

#' @keywords internal
#' @noRd
.check_coord_object <- function(x,
                                check_seqnames = FALSE,
                                check_strand = FALSE) {
  if (!is.data.frame(x)) {
    stop(
      "object must be a data.frame. ",
      "GRanges objects are currently not supported and must be converted ",
      "using e.g. as.data.frame()"
    )
  }
  
  if (!all(c("start", "end") %in% colnames(x))) {
    stop("object must have the columns 'start' and 'end'")
  }
  
  if (check_seqnames) {
    if (!("seqnames" %in% colnames(x))) {
      stop("object must have the column 'seqnames'")
    }
  }
  
  if (check_strand) {
    if (!("strand" %in% colnames(x))) {
      stop("object must have the column 'strand'")
    }
  }
}

#' @keywords internal
#' @noRd
.check_group_var <- function(x, group_var) {
  if (!is.null(group_var)) {
    if (!all(group_var %in% colnames(x))) {
      stop(
        "group_var ('", group_var, "') ",
        "must be a column in object"
      )
    }
  }
}