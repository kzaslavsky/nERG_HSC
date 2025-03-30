#Graphing functions for nERG study
#Date: Feb 25 2024
#Ver: 1.0

library(RColorBrewer)

#set color
color_theme <- brewer.pal(8, "Accent")[c(1,5,6,2,3,4,7,8)]

exp_theme_pt <- function (gp,factorVar = "uni0_or_multimodal1", scale = 1, alfa = 0.5, sz = 3, strk = 1, wd = 0.15 )
{
  gp + 
    geom_point(alpha=alfa, size = sz, stroke = strk, aes_string(fill = factorVar), color = "black", shape = 21, 
               position = position_jitter(width=wd)) +
    #stat_summary(fun.y = mean, geom = "point", alpha = 0.9, size = 5, shape = 16) +
    stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                 geom = "crossbar", color = "black",width = 0.3) +
    stat_summary(fun.data = mean_se, geom = "errorbar", alpha=0.9, lwd=1.25, width = 0.1, color = "black")+
    
    theme(
      axis.title.y = element_text(face = "bold", size = 16*scale, margin = margin(0,10,0,10)),
      axis.title.x = element_text(face = "bold", size = 16*scale, margin = margin(30,0,10,0)),
      axis.text.x = element_text(size = 14*scale, face = "bold", color = "black", hjust = 0.5),
      axis.text.y = element_text(size = 14*scale, face = "bold", color = "black"),
      axis.ticks.x = element_line(color = "black", size = 1),
      strip.text = element_text(size = 16*scale, face = "bold"),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "white"),
      axis.line = element_line(color = "black", size = 1),
      legend.background = element_rect(fill=NA, color = "black"),
      legend.text = element_text(colour="black", size = 14*scale, face = "bold"),
      legend.title = element_text(colour="black", size = 16*scale, face = "bold"),
      legend.position = "right",
      axis.line.y = element_line(color = "black", size = 1),
      axis.line.x = element_line(color = "black", size = 1),
      plot.title = element_text(lineheight=.8, face="bold", size = 20*scale, hjust = 0.5),
      axis.ticks = element_line(color = "black", size = 1)
    )
}

exp_theme <- function (gp,scale = 1, alfa = 0.5 )
{
  gp + 
    theme(
      axis.title.y = element_text(face = "bold", size = 16*scale, margin = margin(0,10,0,10)),
      axis.title.x = element_text(face = "bold", size = 16*scale, margin = margin(30,0,10,0)),
      axis.text.x = element_text(size = 14*scale, face = "bold", color = "black", hjust = 0.5),
      axis.text.y = element_text(size = 14*scale, face = "bold", color = "black"),
      axis.ticks.x = element_line(color = "black", size = 1),
      strip.text = element_text(size = 16*scale, face = "bold"),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "white"),
      axis.line = element_line(color = "black", size = 1),
      legend.background = element_rect(fill=NA, color = "black"),
      legend.text = element_text(colour="black", size = 14*scale, face = "bold"),
      legend.title = element_text(colour="black", size = 16*scale, face = "bold"),
      legend.position = "right",
      axis.line.y = element_line(color = "black", size = 1),
      axis.line.x = element_line(color = "black", size = 1),
      plot.title = element_text(lineheight=.8, face="bold", size = 20*scale, hjust = 0.5),
      axis.ticks = element_line(color = "black", size = 1)
    )
}

exp_theme_hist <- function (gp,scale = 1, alfa = 0.5 )
{
  gp + 
    theme_blank()+
    theme(
      axis.title.y = element_text(face = "bold", size = 16*scale, margin = margin(0,10,0,10)),
      axis.title.x = element_text(face = "bold", size = 16*scale, margin = margin(30,0,10,0)),
      axis.text.x = element_text(size = 14*scale, face = "bold", color = "black", hjust = 0.5),
      axis.text.y = element_text(size = 14*scale, face = "bold", color = "black"),
      axis.ticks.x = element_line(color = "black", size = 1),
      strip.text = element_text(size = 16*scale, face = "bold"),
      panel.border = element_rect(color = "black", size = 2),
      panel.background = element_rect(fill = "white"),
      axis.line = element_line(color = "black", size = 1),
      legend.background = element_rect(fill=NA, color = "black"),
      legend.text = element_text(colour="black", size = 14*scale, face = "bold"),
      legend.title = element_text(colour="black", size = 16*scale, face = "bold"),
      legend.position = "right",
    #  axis.line.y = element_line(color = "black", size = 1),
     # axis.line.x = element_line(color = "black", size = 1),
      plot.title = element_text(lineheight=.8, face="bold", size = 20*scale, hjust = 0.5),
      axis.ticks = element_line(color = "black", size = 1)
    )
}


exp_theme_scatter <- function (gp, factorVar = "Treatment", scale = 1 , alfa = 0.5, shp =21, sz = 5, strk = 1)
{
  gp + 
    geom_point(alpha=alfa, size = sz, stroke = strk, aes_string(fill = factorVar), shape = shp) +
    theme(
      axis.title.y = element_text(face = "bold", size = 16*scale, margin = margin(0,10,0,10)),
      axis.title.x = element_text(face = "bold", size = 16*scale, margin = margin(30,0,10,0)),
      axis.text.x = element_text(size = 14*scale, face = "bold", color = "black", hjust = 0.5),
      axis.text.y = element_text(size = 14*scale, face = "bold", color = "black"),
      axis.ticks.x = element_line(color = "black", size = 1),
      strip.text = element_text(size = 16*scale, face = "bold"),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "white"),
      axis.line = element_line(color = "black", size = 1),
      legend.background = element_rect(fill=NA, color = "black"),
      legend.text = element_text(colour="black", size = 14*scale, face = "bold"),
      legend.title = element_text(colour="black", size = 16*scale, face = "bold"),
      legend.position = "right",
      axis.line.y = element_line(color = "black", size = 1),
      axis.line.x = element_line(color = "black", size = 1),
      plot.title = element_text(lineheight=.8, face="bold", size = 20*scale, hjust = 0.5),
      axis.ticks = element_line(color = "black", size = 1)
    )
}

exp_theme_color <- function (gp)
{
  gp +
    scale_color_manual(values = color_theme[c(3,2)])+
    scale_fill_manual(values = color_theme[c(3,2)])
}



delete <- function(DT, del.idxs) {           # pls note 'del.idxs' vs. 'keep.idxs'
  keep.idxs <- setdiff(DT[, .I], del.idxs);  # select row indexes to keep
  cols = names(DT);
  DT.subset <- data.table(DT[[1]][keep.idxs]); # this is the subsetted table
  setnames(DT.subset, cols[1]);
  for (col in cols[2:length(cols)]) {
    DT.subset[, (col) := DT[[col]][keep.idxs]];
    DT[, (col) := NULL];  # delete
  }
  return(DT.subset);
}

lm_eqn <- function(df){
  m <- lm(df[,1] ~ df[,2], df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}


#function that cna be used for MVN pdf
mvnpdf <- function (x,XBAR,S)
{
  xtemp = matrix(x)
  V = length(XBAR)
  outer_denom = ((2*pi)^V)*(det(S)^.5)
  exp_numerator = t(xtemp-XBAR) %*% solve(S) %*%(xtemp-XBAR)
  mvnpdf = (1/outer_denom)*exp(-1*exp_numerator/2)
}

labs <- c("First-Line", "Second-Line", "DMF")

lmeTable = function(x) {
  fromLme = x
  thetable = summary(fromLme)$tTable
  toRep = dim(thetable)[2]-1
  if(length(fromLme$modelStruct$reStruct)==1) {
    newnames = c(rownames(thetable), c('$\\sigma$', '$\\tau$'))
    thetable = rbind(thetable, 
                     c(sqrt(as.matrix(fromLme$modelStruct$reStruct[[1]]))* fromLme$sigma, rep(NA, toRep)),
                     c(fromLme$sigma, rep(NA, toRep)) 
    )
  } else {
    newnames = c(rownames(thetable), c(names(fromLme$modelStruct$reStruct), '$\\tau$'))
    thetable = rbind(thetable, 
                     cbind(sqrt(diag(as.matrix(fromLme$modelStruct$reStruct[[1]])))* fromLme$sigma, matrix(NA, length(fromLme$modelStruct$reStruct), toRep)),
                     c(fromLme$sigma, rep(NA, toRep)) )
  }
  rownames(thetable) = newnames
  colnames(thetable)[colnames(thetable)=='Value'] = 'MLE'
  
  # if model is ar1
  if('corStruct' %in% names(summary(fromLme$modelStruct))){
    rangeNugget = c(
      summary(fromLme$modelStruct)$corStruct[[1]],
      summary(fromLme$modelStruct)$corStruct[[2]])
    rangeNugget = exp(rangeNugget)
    thetable = rbind(thetable, range=NA, sigmav = NA)
    tauRow = grep("tau\\$$", rownames(thetable))
    sigmaRow = grep("sigma\\$$", rownames(thetable))
    thetable['range','MLE'] = rangeNugget[1]
    tausq = rangeNugget[2] * thetable[tauRow, 'MLE']^2
    thetable['sigmav','MLE'] = sqrt( 
      thetable[tauRow, 'MLE']^2 - tausq)	
    thetable[tauRow, 'MLE'] =  sqrt(tausq)
    rownames(thetable) = gsub("sigmav", "$\\\\sigma_V$",
                              rownames(thetable))
    rownames(thetable)[sigmaRow] = '$\\sigma_U$'
  }
  thetable
}


lm_eqn <- function(df){
  m <- lm(y ~ x, df)
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(as.numeric(m$coefficients[1]),digits = 2),
                        b = format(as.numeric(m$coefficients[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 2)))
  as.character(as.expression(eq));     
}


# hgvsParseR is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# hgvsParseR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with hgvsParseR.  If not, see <https://www.gnu.org/licenses/>.

#' HGVS Parser 
#'
#' Parses HGVS strings
#' @param strings A character vector containing the HGVS strings
#' @param aacode allowed values: 1, 3, or NA. Determines whether 1-letter codes or 3-letter codes should be forced. NA uses input format.
#' @return A \code{data.frame} with the following columns: 
#' @keywords HGVS parsing
#' @export
#' @examples
#' result <- parseHGVS(c("g.1318G>T","c.123_125inv","p.R123_L152del"))



strReverse <- function(x)
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")


