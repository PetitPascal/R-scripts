#-------------------------------------------------------------------------------
## Reproducible and generalizable function enabling to save/export any ggplot2 object or a list of ggplots with automatic size
#-------------------------------------------------------------------------------

#----------------------------------------------------------------
#### Configurations ####

## General R environment step-up

rm(list=ls()) # clearing the environment
gctorture(FALSE) # disabling memory torture
options(stringsAsFactors=FALSE) # no automatic conversion of characters into factors

## R packages

required_pkgs<-c("tidyverse","conflicted","patchwork")

is_installed<-required_pkgs %in% rownames(installed.packages(all.available=TRUE))
if(any(is_installed==FALSE)){
  install.packages(required_pkgs[!is_installed],repos="http://cran.us.r-project.org")
}
invisible(lapply(required_pkgs, library, character.only=TRUE))

## Preventing package conflicts

conflict_prefer("select","dplyr")
conflict_prefer("filter","dplyr")
conflict_prefer("slice","dplyr")

## Setting the working directory

here::here("Figure export")

#----------------------------------------------------------------
#### Creating the function ####

#- - - - - -
## Saving figures (any ggplot2 object, or a list of ggplots) with automatic size
#
# - ext: one or several formats, e.g. c("png","pdf","tiff"); one file per format (name_<date>.<ext>)
# - compression: TIFF compression (e.g., "lzw"); only used for tiff files
# - discrete axes: the size is computed from the number and length of the labels and from the text written inside the tiles (geom_text/geom_label)
# - continuous axes: panel_width x panel_height (cm)
# - facets with free space keep their relative sizes; fixed-ratio coordinates (coord_fixed, coord_sf) keep their aspect ratio
# - one plot: the panel size is set directly in the plot layout and the final size is measured on this layout (panels + axes + strips + legend), widened if the title/subtitle/caption is wider
# - several plots: one single legend, placed at legend_position ("top" = above all plots)
# - a patchwork object cannot be measured: width/height are used (default 25 x 18 cm)
# - width, height: optional (cm); override the automatic size
# - ...: other arguments passed to ggplot2::ggsave

save_plot<-function(plot, name, ncol=1, dpi=300, ext="png", panel_width=9, panel_height=7, legend_position=NULL, margin=0.5,
                    width=NULL, height=NULL, compression="lzw", ...){
  plots<-if(inherits(plot, "ggplot")) list(plot) else plot
  ncol<-min(ncol, length(plots))
  nrow_p<-ceiling(length(plots)/ncol)
  ext<-unique(tolower(ext))
  files<-paste0(name, "_", Sys.Date(), ".", ext) # one file per format
  
  # null device: needed by grid to measure text, without creating a file
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add=TRUE)
  
  pt_cm<-0.03528 # 1 point in cm
  
  # width (cm) of the widest text of a vector
  text_cm<-function(x, fontsize){
    x<-as.character(x[!is.na(x)])
    if(length(x)==0) return(0)
    grid::convertWidth(grid::grobWidth(grid::textGrob(x, gp=grid::gpar(fontsize=fontsize))), "cm", valueOnly=TRUE)
  }
  # value of a theme element of the plot (e.g., font size of the axis text)
  theme_el<-function(p, el, field, default){
    th<-tryCatch(ggplot2::theme_get() + p$theme, error=function(e) ggplot2::theme_get())
    v<-tryCatch(ggplot2::calc_element(el, th)[[field]], error=function(e) NULL)
    if(is.null(v) || length(v)==0 || is.na(v[1])) default else v
  }
  
  # panel area of one plot: returns its size (cm) and its layout (gtable) with the panels set to this size
  panel_dims<-function(p){
    b<-ggplot2::ggplot_build(p)
    g<-ggplot2::ggplotGrob(p)
    lay<-b$layout$layout
    pan<-g$layout[grepl("^panel", g$layout$name), ]
    cols<-sort(unique(pan$l)) # columns and rows of the layout occupied by panels
    rows<-sort(unique(pan$t))
    
    # with coord_flip, the x scale is drawn vertically
    flipped<-inherits(p$coordinates, "CoordFlip")
    h_ax<-if(flipped) "y" else "x"
    v_ax<-if(flipped) "x" else "y"
    sc_h<-b$layout[[paste0("panel_scales_", h_ax)]]
    sc_v<-b$layout[[paste0("panel_scales_", v_ax)]]
    idx_h<-lay[[paste0("SCALE_", toupper(h_ax))]]
    idx_v<-lay[[paste0("SCALE_", toupper(v_ax))]]
    
    # text layers (all labels need room, including repelled ones)
    cls<-sapply(p$layers, function(l) class(l$geom)[1])
    is_txt<-grepl("Text|Label", cls)
    txt<-if(any(is_txt)) as.character(unlist(lapply(b$data[which(is_txt)], function(d) d$label))) else character(0)
    has_txt<-length(txt)>0
    
    # text written inside the tiles (plain geom_text / geom_label, not repelled): the slots of discrete axes must be large enough to contain it
    is_tile<-cls %in% c("GeomText","GeomLabel")
    tile_txt<-if(any(is_tile)) as.character(unlist(lapply(b$data[which(is_tile)], function(d) d$label))) else character(0)
    tile_txt<-tile_txt[!is.na(tile_txt) & tile_txt!=""]
    tile_fs<-if(any(is_tile)) max(unlist(lapply(b$data[which(is_tile)], function(d) if("size" %in% names(d)) d$size else 3.88)), na.rm=TRUE) * ggplot2::.pt else 11
    tile_w<-if(length(tile_txt)>0) 1.15 * text_cm(unlist(strsplit(tile_txt, "\n", fixed=TRUE)), tile_fs) + 0.6 else 0 # widest line (cm): +15% and padding
    tile_h<-if(length(tile_txt)>0) (max(stringr::str_count(tile_txt, "\n")) + 1) * tile_fs * pt_cm * 1.4 + 0.5 else 0 # tallest text (cm)
    
    fs_x<-theme_el(p, "axis.text.x", "size", 11)
    fs_y<-theme_el(p, "axis.text.y", "size", 11)
    ang_x<-theme_el(p, "axis.text.x", "angle", 0)
    lab_of<-function(s) tryCatch(as.character(s$get_labels()), error=function(e) as.character(s$get_limits()))
    
    # horizontal size of one panel: discrete = one slot per label (widest of: label, text inside the tile); continuous = panel_width
    size_h<-function(s){
      if(s$is_discrete()){
        lab<-lab_of(s)
        fs<-if(flipped) fs_y else fs_x
        per<-if(!flipped && ang_x>=45) max(0.9, fs * pt_cm * 1.8) else max(0.9, text_cm(lab, fs) + 0.5)
        max(4, length(lab) * max(per, tile_w))
      }else{
        max(panel_width, if(has_txt) 3 * text_cm(txt, 11) else 0)
      }
    }
    # vertical size of one panel: discrete = one slot per label (height of its lines, or of the text inside the tile); continuous = panel_height
    size_v<-function(s){
      if(s$is_discrete()){
        lab<-lab_of(s)
        fs<-if(flipped) fs_x else fs_y
        n_l<-stringr::str_count(lab, "\n") + 1
        per<-pmax(if(has_txt) 1.1 else 0.8, n_l * fs * pt_cm * 1.25 + 0.35, tile_h)
        max(3, sum(per))
      }else{
        panel_height
      }
    }
    w<-sapply(seq_len(max(lay$COL)), function(i) max(sapply(unique(idx_h[lay$COL==i]), function(k) size_h(sc_h[[k]]))))
    h<-sapply(seq_len(max(lay$ROW)), function(i) max(sapply(unique(idx_v[lay$ROW==i]), function(k) size_v(sc_v[[k]]))))
    if(length(cols)!=length(w)) w<-rep(max(w), length(cols))
    if(length(rows)!=length(h)) h<-rep(max(h), length(rows))
    
    # facets with free space (unequal null units) keep their relative sizes, otherwise all panels are equal
    keep_prop<-function(sizes, units){
      nul<-tryCatch(as.numeric(units), error=function(e) rep(1, length(units)))
      if(length(nul)>1 && diff(range(nul))>1e-8) nul/max(nul) * max(sizes) else rep(max(sizes), length(sizes))
    }
    w<-keep_prop(w, g$widths[cols])
    h<-keep_prop(h, g$heights[rows])
    
    # fixed aspect ratio (coord_fixed, coord_sf...): the height follows the width
    free<-tryCatch(p$coordinates$is_free(), error=function(e) TRUE)
    if(!free){
      asp<-tryCatch(sum(as.numeric(g$heights[rows]))/sum(as.numeric(g$widths[cols])), error=function(e) 1)
      w<-rep(panel_width/length(w), length(w))
      h<-rep(sum(w) * asp/length(h), length(h))
      if(sum(h)>2 * panel_height){
        h<-h * 2 * panel_height/sum(h)
        w<-w * (sum(h)/asp)/sum(w)
      }
    }
    
    # panels set to the computed size in the layout
    g$widths[cols]<-grid::unit(w, "cm")
    g$heights[rows]<-grid::unit(h, "cm")
    list(w=sum(grid::convertWidth(g$widths[min(cols):max(cols)], "cm", valueOnly=TRUE)),
         h=sum(grid::convertHeight(g$heights[min(rows):max(rows)], "cm", valueOnly=TRUE)),
         grob=g)
  }
  
  if(any(sapply(plots, inherits, what="patchwork"))){
    # already assembled patchwork: cannot be measured
    message("A patchwork object cannot be measured: width/height used (default 25 x 18 cm); use width= and height= (cm) to change them")
    fig<-if(length(plots)==1) plots[[1]] else patchwork::wrap_plots(plots, ncol=ncol, guides="collect")
    if(!is.null(legend_position)) fig<-fig & ggplot2::theme(legend.position=legend_position)
    if(is.null(width)) width<-25
    if(is.null(height)) height<-18
    
  }else if(length(plots)==1){
    # one plot: panels set directly in the layout of the plot, size measured on this layout (panels + axes + strips + legend)
    p<-plots[[1]]
    if(!is.null(legend_position)) p<-p + ggplot2::theme(legend.position=legend_position)
    d<-panel_dims(p)
    message("Panel area computed: ", round(d$w, 1), " cm wide x ", round(d$h, 1), " cm high")
    fig<-d$grob
    w_layout<-sum(grid::convertWidth(fig$widths, "cm", valueOnly=TRUE))
    h_layout<-sum(grid::convertHeight(fig$heights, "cm", valueOnly=TRUE))
    # title, subtitle and caption are not part of the column widths: the figure must be at least as wide as the longest one
    w_text<-max(c(0, sapply(c("title","subtitle","caption"), function(l){
      s<-tryCatch(p$labels[[l]], error=function(e) NULL)
      if(is.null(s) || length(s)==0 || is.na(s[1])) return(0)
      1.1 * text_cm(strsplit(as.character(s), "\n", fixed=TRUE)[[1]], theme_el(p, paste0("plot.", l), "size", 11))
    })))
    if(is.null(width)) width<-max(w_layout * 1.02, w_text) + margin
    if(is.null(height)) height<-h_layout * 1.02 + margin
    
  }else{
    # several plots: panel area of each column and each row of plots = largest of its plots, assembled with patchwork
    sizes<-lapply(plots, panel_dims)
    idx<-matrix(seq_len(nrow_p * ncol), nrow=nrow_p, ncol=ncol, byrow=TRUE)
    w_col<-sapply(seq_len(ncol), function(j) max(sapply(idx[, j][idx[, j]<=length(plots)], function(k) sizes[[k]]$w)))
    h_row<-sapply(seq_len(nrow_p), function(i) max(sapply(idx[i, ][idx[i, ]<=length(plots)], function(k) sizes[[k]]$h)))
    message("Panel area computed: ", paste0(round(w_col, 1), collapse=" / "), " cm wide x ", paste0(round(h_row, 1), collapse=" / "), " cm high (per column / per row)")
    # one single legend: kept on the first plot only (all plots of a figure must share the same color scale), removed from the others
    if(!is.null(legend_position)){
      plots<-lapply(seq_along(plots), function(k) plots[[k]] + ggplot2::theme(legend.position=if(k==1) legend_position else "none"))
    }
    # legend on top of everything (above the titles of the first row): a dedicated first row is reserved for the collected legend
    top_legend<-identical(legend_position, "top")
    if(top_legend){
      # height of the legend (cm), measured on the first plot
      leg_h<-tryCatch({
        gl<-ggplot2::ggplotGrob(plots[[1]] + ggplot2::theme(legend.position="top"))
        top_panel<-min(gl$layout$t[grepl("^panel", gl$layout$name)])
        i<-which(grepl("guide-box", gl$layout$name) & gl$layout$t<top_panel)
        if(length(i)==0) 1.5 else sum(grid::convertHeight(gl$heights[unique(gl$layout$t[i])], "cm", valueOnly=TRUE))
      }, error=function(e) 1.5)
      # row 1 = area of the collected legend (all columns), then one area per plot (filled row by row)
      areas<-Reduce(c, c(list(patchwork::area(1, 1, 1, ncol)),
                         lapply(seq_along(plots), function(k) patchwork::area((k - 1) %/% ncol + 2, (k - 1) %% ncol + 1))))
      fig<-patchwork::wrap_plots(c(list(patchwork::guide_area()), plots), design=areas, widths=grid::unit(w_col, "cm"),
                                 heights=grid::unit(c(leg_h, h_row), "cm"), guides="collect")
    }else{
      fig<-patchwork::wrap_plots(plots, ncol=ncol, widths=grid::unit(w_col, "cm"), heights=grid::unit(h_row, "cm"), guides="collect")
    }
    # total size = size of the assembled figure (+2% and a margin)
    gt<-patchwork::patchworkGrob(fig)
    if(is.null(width)) width<-sum(grid::convertWidth(gt$widths, "cm", valueOnly=TRUE)) * 1.02 + margin
    if(is.null(height)) height<-sum(grid::convertHeight(gt$heights, "cm", valueOnly=TRUE)) * 1.02 + margin
  }
  
  # one file per format; compression is only passed to tiff devices
  for(k in seq_along(files)){
    args<-c(list(plot=fig, filename=files[k], width=width, height=height, units="cm", dpi=dpi, bg="white", limitsize=FALSE), list(...))
    if(ext[k] %in% c("tiff","tif")){
      args$compression<-compression
      args$device<-if(requireNamespace("ragg", quietly=TRUE)) ragg::agg_tiff else "tiff"
    }
    do.call(ggplot2::ggsave, args)
    message("Figure saved: ", files[k], " (", round(width, 1), " x ", round(height, 1), " cm)")
  }
  invisible(c(width=width, height=height))
}

#----------------------------------------------------------------
#### Example of use ####

#- - - - - 
## Creating plots

plot1<-ggplot(airquality,aes(x=factor(Day),y=Ozone))+
  geom_boxplot()+
  scale_x_discrete("Day")+
  scale_y_continuous("Ozone concentration")+
  theme_bw()

plot2<-ggplot(airquality, aes(Day, Ozone)) +
  geom_smooth(stat = 'summary', 
              fun.data = function(y) data.frame(ymin = quantile(y, .05), y = mean(y), ymax = quantile(y, .95)))+
  theme_bw()

plot3<-ggplot(airquality, aes(x=Month, y=Day, fill=Ozone)) + 
  geom_tile(color="black")+ 
  scale_fill_distiller(palette = "RdBu")+
  theme_bw()

plot4<-airquality %>% mutate(Wind=case_when(Wind<10~"Wind<10",
                                            Wind>15~"Wind>15",
                                            T~"Wind: 10-15")) %>%
  ggplot() + geom_boxplot(aes(x=factor(Month),y=Ozone)) +
  scale_y_continuous("O3 concentration",breaks=seq(0,200,25), labels=seq(0,200,25), limits=c(0,200)) +
  scale_x_discrete("Mois",breaks=c(5:9), labels=c("May","Jun.","Jul.","Aug.","Sept."))+
  theme_bw()+
  facet_grid(.~Wind)

plot5<-plot2+plot3/plot4

#- - - - - 
## Exporting the graph

# export in png
save_plot(plot1, "Plot 1", ncol=3, legend_position="top", ext="png")

# export in pdf
save_plot(plot2, "Plot 2", ncol=3, legend_position="top", ext="pdf")

# export in tiff
save_plot(plot3, "Plot 3", ncol=3, legend_position="top", ext="tiff", dpi=600, compression="lzw")

# export in several formats simultaneously
save_plot(plot4, "Plot 4", ncol=3, legend_position="top", ext=c("png","pdf","tiff"))
save_plot(plot5, "Plot 5", ncol=3, legend_position="top", ext=c("png","pdf","tiff"))
