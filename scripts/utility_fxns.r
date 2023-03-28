
ccdf <- function(x)
{
	1 - ecdf(sort(x))(sort(x))
}

#########
image_table <- function(x, title=character(0), col=c(), zero_col=NULL, colrange=c(), show_xlab=T, show_ylab=T, xlabels=c(), ylabels=c(), xlabels_at=c(), ylabels_at=c(), grid_at=c(), grid_col="black", cex=1, nshades=1000, square=F, colorkey=T, cor=F)
{
	if (show_xlab) {
		if (is.null(xlabels)) {
			xlabels <- colnames(x)
		}

		if (is.null(xlabels_at)) {
			xlabels_at <- c(1:length(xlabels))
		}
	}
	if (show_ylab) {
		if (is.null(ylabels)) {
			ylabels <- rownames(x)[nrow(x):1]
		}

		if (is.null(ylabels_at)) {
			ylabels_at <- c(1:length(ylabels))
		}
	}

	# color stuff
	# sequential favorites
	pubu <- c("#fff7fb","#ece7f2","#d0d1e6","#a6bddb","#74a9cf","#3690c0","#0570b0","#045a8d","#023858")
	pubugn <- c("#fff7fb","#ece2f0","#d0d1e6","#a6bddb","#67a9cf","#3690c0","#02818a","#016c59","#014636")
	bupu <- c("#f7fcfd","#e0ecf4","#bfd3e6","#9ebcda","#8c96c6","#8c6bb1","#88419d","#810f7c","#4d004b")
	gnbu <- c("#f7fcf0","#e0f3db","#ccebc5","#a8ddb5","#7bccc4","#4eb3d3","#2b8cbe","#0868ac","#084081")

	ylgnbu <- c("#edf8b1", "#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#253494", "#081d58")
	reds <- c("#fff5f0","#fee0d2","#fcbba1","#fc9272","#fb6a4a","#ef3b2c","#cb181d","#a50f15","#67000d")
	oranges <- c("#fff5eb","#fee6ce","#fdd0a2","#fdae6b","#fd8d3c","#f16913","#d94801","#a63603","#7f2704")
	ylorbr <- c("#ffffe5","#fff7bc","#fee391","#fec44f","#fe9929","#ec7014","#cc4c02","#993404","#662506")

	# diverging favorites
	rdbu <- c("#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac")
	rdylgn <- c("#a50026","#d73027","#f46d43","#fdae61","#fee08b","#ffffbf","#d9ef8b","#a6d96a","#66bd63","#1a9850","#006837")

	# set the color ramp, and reverse if needed
	if (is.null(col)) {
		if (!is.null(zero_col)) {
			pubu <- pubu[2:length(pubu)]
		}
		color.ramp <- colorRampPalette(pubu)(nshades)
	}
	else {
		colors = get(sub("\\*","",col))
		if (!is.null(zero_col)) {
			colors <- colors[2:length(colors)]
		}
		color.ramp <- colorRampPalette(colors)(nshades)
		if (grepl("\\*", col)) {
			color.ramp <- rev(color.ramp)
		}
	}
	
	if (!is.null(zero_col)) {
		color.ramp[1] <- zero_col
	}

	# if no color range specified, choose according to correlation or just min/max
	if (is.null(colrange)) {
		if (cor==T) {
			bound <- max(abs(x[!is.na(x) & x!=1]))
			color.levels <- seq(-bound, bound, length=nshades)
			colrange <- c(-bound, bound)
		}
		else {
			color.levels <- seq(min(x[!is.na(x)]), max(x[!is.na(x)]), length=nshades)
			colrange <- c(color.levels[1], color.levels[length(color.levels)]) 
		}
	}
	else {
		color.levels <- seq(colrange[1], colrange[2], length=nshades)
	}

	# plot and colorkey dimensions
	if (colorkey) {
		if (square) {
			map.rel.width <- 8
			pty.hm="m"
		}
		else {
			map.rel.width <- 18
			pty.hm="m"
		}
		pty.ck="m"
		layout(matrix(data=c(rep(1,map.rel.width),rep(2,2)), nrow=1, ncol=map.rel.width+2), heights=c(1,1))
	}
	else {
		pty.hm="m"
	}

	if (length(title)==0) {
		mar.heatmap <- c(2,5*cex,4*cex,0)
		mar.colorkey <- c(2,1,4*cex,4*cex)
	}
	else {
		mar.heatmap <- c(2,5*cex,10*cex,0)
		mar.colorkey <- c(2,1,10*cex,4*cex)
	}

	# matrix heatmap
	par(mar=mar.heatmap, pty=pty.hm)
	if (is.null(colrange)) {
		image(x=1:ncol(x), y=1:nrow(x),
			# z=t(x)[nrow(x):1,nrow(x):1],
			z=t(x),
			ylim=c(nrow(x),1),
			xlab="", ylab="", axes=F,
			col=color.ramp)
	}
	else {
		image(x=1:ncol(x), y=1:nrow(x),
			# z=t(x)[nrow(x):1,nrow(x):1],
			z=t(x),
			ylim=c(nrow(x),1),
			xlab="", ylab="", axes=F,
			col=color.ramp, zlim=colrange)
	}
	if (!is.null(grid_at)) {
		abline(h=grid_at, v=grid_at, col=grid_col)
	}
	box()
	if (show_xlab) {
		axis(side=3, at=xlabels_at, labels=xlabels, cex.axis=2*cex)
	}
	if (show_ylab) {
		axis(side=2, at=ylabels_at, labels=ylabels, las=HORIZONTAL<-1, cex.axis=2*cex)
	}

	if (length(title)!=0) {
		title(title, cex.main=cex*2.5)
	}

	# color key
	if (colorkey) {
		par(mar=mar.colorkey, pty=pty.ck)
		# par(pty=pty.ck)
		image(1, color.levels,
			matrix(data=color.levels, ncol=length(color.levels), nrow=1),
			col=color.ramp, xlab="",ylab="", axes=F)
		box()
		axis(4, cex.axis=cex*2, las=2)
	}
}

#########
# from stack overflow

simple_cap <- function(x)
{
	# s <- strsplit(x, " ")[[1]]
	paste(toupper(substring(x,1,1)), substring(x,2), sep="", collapse=" ")
}

#########

# get_lower_tri_idx <- function(x, incl_diag=T)
# {
# 	nrows <- nrow(x)
# 	z <- seq(nrows)

# 	if (incl_diag) {
# 		row <- unlist(lapply(1:nrows, function(x) x:nrows), use.names=F)
# 		col <- rep(z[-length(z)], times=rev(tail(z,-1)))
# 	}
# 	else {
# 		row <- unlist(lapply(2:nrows, function(x) x:nrows), use.names=F)
# 		col <- rep(z[-length(z)], times=rev(tail(z,-1))-1)
# 	}

# 	cbind(row, col)
# }

########
# Wrapper to AT K-means clustering (impl by C) that supports NA values
TGLKMeans_wrapper=function(data, fn, k)
{
  write.table(data, fn, sep="\t", quote=F, col.names=NA)
  system(sprintf("/Users/Arnau/Dropbox/Cris_RNAseq/TGLKMeans_static %s %s euclid -allow_nas=1 &> %s.log", fn, k, fn))
  km = list()
  m.k = read.table(paste(fn, "kclust", sep="."), header=T)
  m.k$clust = m.k$clust + 1
  m.c = read.delim(paste(fn, "center", sep="."), header=F)
  km$size = tapply(m.k$id, m.k$clust, length)
  km$cluster = m.k$clust
  names(km$cluster) = m.k$id
  km$centers = m.c[, 2:ncol(m.c)]
  km
} 

###################
# Input: list of cells, chromosome name, quantile to filter counts on, and which contacts to count (options: self, X2M, far and trans)
# Output: List of filtered cells and the minimum number of contacts they have on the chromosome
# Note that the contacts are doubled, so x-y and y-x are considered 2 contacts
filter_by_contacts_per_chrom <<- function(cells, chrom, quant, filter_by=c('X2M', 'far'))
{
  chrom = ifelse(length(grep("^chr", chrom)) > 0, chrom, paste("chr", chrom, sep=""))
  cols = paste(filter_by, chrom, sep="_")
  x = rowSums(sch_chrom_stat[cells, cols])
  n = quantile(x, quant) - 1
  x = x[x >= n]
  list(cells=names(x), min_contacts=n)
} 

####################
# data with k columns, plot k*k matrix with scatter plots between columns and correlation 
plot_scatter_cols <<- function(data, name, width, height) {
  png(name, width=width, height=height)
  layout(matrix(1:ncol(data)**2, ncol(data), ncol(data)))
  cols = colorRampPalette(nettas_trio)(201)
  par(mar=c(2,2,1,1))
  for (j in 1:ncol(data)) {
    for (i in 1:ncol(data)) {
      pc = round(cor(data[,i], data[,j], use="pair"), 2)
      
      if (i == j) {
        plot.new()
        plot.window(0:1, 0:1)
        rect(0, 0, 1, 1, col='gray')
        text(x=0.1, y=0.75, labels=colnames(data)[i], pos=4, cex=4)
        text(x=0.1, y=0.25, labels=signif(var(data[,i]), 2), pos=4, cex=3)
      }
      if (i > j) {
        plot(data[,i], data[,j], pch=19, cex=0.7, xlab="", ylab="", main="") # xlim=range(data), ylim=range(data))
        grid(col='gray')
      }
      if (j > i) {
        plot.new()
        plot.window(0:1, 0:1)
        rect(0, 0, 1, 1, col=cols[101 + 100 * pc])
        text(x=0.5, y=0.5, labels=pc, cex=3, col='white')
        
      }
    }
  }
  dev.off()
}

####################
# Plot scatter of x vs y (vectors) and color by columns of z (ncol(z) subplots)
plot_colored_scatters <<- function(x, y, z, z_discrete, ofn, xlab=NULL, ylab=NULL)
{
  png(ofn, width=ncol(z)*250, height=260)
  layout(t(1:ncol(z)))

  for (i in 1:ncol(z)) {
    v = z[,i]
    ind = !is.na(x) & !is.na(y) & !is.na(v)
    
    if (z_discrete[i])  {
      cols = qualitative_colors[v]
    }
    else {
      breaks = c(min(v, na.rm=T), median(v, na.rm=T), max(v, na.rm=T))
      cols = colorRampPalette(c("blue", "black", "yellow"))(256)[vals_to_cols(v, breaks, 128)]
    }

    plot(x[ind], y[ind], pch=19, cex=1.25, main=colnames(z)[i], xlab=xlab, ylab=ylab, col=cols[ind])
  }
  dev.off()
}

###################
vals_to_cols=function(vals, breaks, ncols=256)
{
  min = breaks[1]
  max = breaks[length(breaks)]
  n = length(breaks)-1
  cols = rep(-1, length(vals))
  for (i in 1:n)
  {
    ind = !is.na(vals) & (breaks[i] <= vals) & (vals <= breaks[i+1])
    if (!any(ind))
      next
    # normalize to [0,1]
    cols[ind] = (vals[ind] - breaks[i]) / (breaks[i+1] - breaks[i])
    # normalize to [i*ncols,i*(ncols+1)]
    cols[ind] = (i-1)*ncols + cols[ind]*(ncols-1) + 1
    # round
    cols[ind] = round(cols[ind])
  }
  cols[cols == -1] = NA
  return (cols)
}

########################
plot_legend=function(ofn, zlim, crp, name)
{
  png(ofn, width=3*length(crp), height=150)
  image(as.matrix(seq(zlim[1], zlim[2], length=length(crp))), zlim=zlim, col=crp, xaxt='n', yaxt='n', main=name)
  axis(1, at=c(0,1), labels=zlim)
  dev.off()
}

########################
n2str = function(n, digits=2)
{
  f    = ifelse(n >= 1e9, 1e9, ifelse(n >= 1e6, 1e6, ifelse(n >= 1e3, 1e3, 1)))
  suff = ifelse(n >= 1e9, "G", ifelse(n >= 1e6, "M", ifelse(n >= 1e3, "K", "")))

  ifelse(digits > 0, sprintf(paste("%.", digits, "f%s", sep=""), n/f, suff), paste(n/f, suff, sep=""))
}

########################## addalpha() from https://github.com/mylesmharrison/colorRampPaletteAlpha/blob/master/colorRampPaletteAlpha.R
addalpha <- function(colors, alpha=1.0) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}

############
load_table_if_exists <- function(fn, header=T, sep="\t") {
  out = NULL
  if (file.exists(fn)) {
    out = read.table(fn, header=header, sep=sep)
  }
  out
}

############
cyclic_rollapply <- function(x, width, by=1, func_str="mean") {
  stopifnot(width %% 2 == 1)
  hsw = (width-1)/2
  if (is.vector(x)) {
    l = length(x)
    y = x[c((l - hsw + 1):l, 1:l, 1:hsw)]
  }
  else {
    l = ncol(x)
    y = t(x[,c((l - hsw + 1):l, 1:l, 1:hsw)])
  }
  
  if (func_str == "mean") {
    y = rollmean(y, k=width, fill=NA)
  }
  else if (func_str == "sum") {
    y = rollsum(y, k=width, fill=NA)
  }
  else if (func_str == "sd") {
    y = rollapply(y, width=width, FUN=sd, partial=T)
  }
  else if (func_str == "length") {
    y = rollapply(y, width=width, FUN=length, partial=T)
  }

  if (is.vector(x)) {
    y = y[seq(hsw + 1, hsw + l,  by=by)]
  }
  else {
    y = t(y[ seq(hsw + 1, hsw + l,  by=by), ])
  }
  y
}

cyclic_rollmean <- function(x, width, by=1) {
  cyclic_rollapply(x, width, by, "mean")
}

cyclic_rollsum <- function(x, width, by=1) {
  cyclic_rollapply(x, width, by, "sum")
}

##############
# intervals manipulation (from Netta)
intervals.centers <- function(inv){
  inv[,2:3]<-floor((inv[,2]+inv[,3])/2)
  inv[,3]<-inv[,3]+1
  return(inv)
}
intervals.2d.centers <- function(inv){
  inv[,2:3]<-floor((inv[,2]+inv[,3])/2)
  inv[,3]<-inv[,3]+1
  inv[,5:6] = floor((inv[,5]+inv[,6])/2)
  inv[,6] = inv[,5] + 1
  return(inv)
}

intervals.expand <- function(inv,expansion=100){
  inv[,2]<-inv[,2]-expansion
  inv[,3]<-inv[,3]+expansion
  return(gintervals.force_range(inv))
}

intervals.2d.expand <- function(inv,expansion1, expansion2){
  inv[,2]<-inv[,2]-expansion1
  inv[,3]<-inv[,3]+expansion1
  inv[,5]<-inv[,5]-expansion2
  inv[,6]<-inv[,6]+expansion2
  return(gintervals.force_range(inv))
}


intervals.size <- function(inv){
  return(sum(inv[,3]-inv[,2]))
}

intervals.normalize <- function(inv, size) {
  centers <- intervals.centers(inv)
  centers$end <- centers$end-1;
  return(intervals.expand(centers, floor(size/2)))
}

###########
.plot_points_and_trend <- function(x, ofn, cols, width_per_panel=500, height_per_panel=110, sw=101, col_by_cond=T, collapse_q=c(0.005, 0.995), cex=0.4, glob_trend=NULL, glob_col=NA, glob_ylim=NULL, lcol="blue", pcol="black", vertical=T, plot_trends=F, pointsize=7, fres=120, trend_cols=mini_qualitative_colors)
{

  np = ncol(x) + ifelse(plot_trends, 1, 0)
  
  if (vertical) {
    png(ofn, width_per_panel, height_per_panel * np, pointsize=pointsize, res=fres)
    par(mar=c(1,4,0,1))
    layout(matrix(1:np, np, 1))
  }
  else {
    png(ofn, width_per_panel * np, height_per_panel, pointsize=pointsize, res=fres)
    par(mar=c(4,2,0,0))
    layout(matrix(1:np, 1, np))
  }

  trends = NULL
  for (i in 1:ncol(x)) {
    v = x[,i]
    
    if (col_by_cond) {
      col = cols
    }
    else {
      col = rep(pcol, length(v))
    }

    if (is.null(glob_ylim)) {
      qs = quantile(v, collapse_q, na.rm=T)
    }
    else {
      qs = glob_ylim
    }

    ind = v < qs[1] | v > qs[2]
    v[v < qs[1]] = qs[1]
    v[v > qs[2]] = qs[2]

    if (is.null(glob_ylim)) {
      ylim = range(v, na.rm=T)
    }
    else {
      ylim = glob_ylim
    }
      
    col[ind] = 'red'

    plot(v, pch=19, cex=cex, main="", ylab="", xaxt='n', xlab="", col=col, ylim=ylim)
    if (vertical) {
      title(ylab=colnames(x)[i])
    }
    else {
      title(xlab=colnames(x)[i])
    }
      
    if (!is.null(glob_trend)) {
      lines(glob_trend, col=glob_col, lwd=2)
    }
    vs = cyclic_rollmean(v, sw)
    lines(vs, col=lcol, lwd=2)
    grid(col='black', nx=NA, ny=NULL)

    trends = rbind(trends, vs)
  }


  if (plot_trends) {
    plot(trends[1,], type='l', ylim=range(trends), col=NA, main="", xlab="", ylab="")
    apply(cbind(trend_cols[1:nrow(trends)], trends), 1, function(y) { lines(y[-1], col=y[1], lwd=2) })
    grid(col='black')
    legend("topleft", legend=colnames(x), col=trend_cols[1:nrow(trends)], lwd=2, bty='n', cex=0.8, ncol=2)
  }
  
  dev.off()  
}

