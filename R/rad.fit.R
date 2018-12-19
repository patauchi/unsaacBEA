#if (getRversion() >= "2.15.1") { utils::globalVariables(c("aes"))}
#' Rango Abundance Distribution Fit
#'
#' \code{Rango Abundance Curve} returns the abundance models.
#'
#' This is a generic function: methods can be defined for it directly
#' or via the \code{\link{Summary}} group generic. For this to work properly,
#' the arguments \code{...} should be unnamed, and dispatch is on the
#' first argument.
#'
#' @param Vdata Two Community with Abundance of species.
#' @param name 
#' @param plot
#' @param method 
#' @param mixing 
#' 
#' @return If all inputs are integer and logical, then the output
#'   will be an integer. If integer overflow
#'   \url{http://en.wikipedia.org/wiki/Integer_overflow} occurs, the output
#'   will be NA with a warning. Otherwise it will be a length-one numeric or
#'   complex vector.
#'
#'   Zero-length vectors have sum 0 by definition. See
#'   \url{http://en.wikipedia.org/wiki/Empty_sum} for more details.
#' @examples
#' a<-c(2,3,4,5,4,5,10)
#' b<-c(0,0,0,1,5,6,50)
#' c<-data.frame(a,b)
#' c<-t(c)
#' beta.diversity(as, indice='jaccard')
#' 
#' @importFrom dplyr %>%
#' @import ggplot2

 rad.fit <- function(Vdata, name=FALSE, plot = F, method='Log10', mixing = FALSE, ...){
   
   # Conditionals
  if(!is.factor(Vdata[[1]]) == TRUE)
    stop('First column is factor class')
   # Extracting dataset
   n <- ncol(Vdata) - 1 # number of communities
   
  radCom <- list()
  #x <- as_tibble(x)
  
  for (i in 1:n) {
    
    k <- i + 1
    mm <- Vdata[,c(1,k)]
      
    mm <- dplyr::filter(mm, mm[[2]] != 0); names(mm) <- c('Species', 'Commun')
    ns <- 1:nrow(mm)
    dat.Ab <- dplyr::mutate(mm, Re.Ab = Commun/sum(Commun), log2pi = log2(Re.Ab), 
                     log10pi = log10(Re.Ab), Community = paste0('Community_', i))
    dat.Ab <- dplyr::arrange(dat.Ab, -Commun)
    dat.Ab <- dplyr::mutate(dat.Ab, Ranks = ns)
    

    radCom[[i]] <- dat.Ab
   
  }
  
  if(plot == FALSE) {
    return(radCom)
  } else{
    
    radCom  <- dplyr::bind_rows(radCom)    
    
    METHODS <- c("Abundance", "Log2","Log10","Rel.Abun")
    
    method <- match.arg(method, METHODS)
    
    if(method == "Abundance") {radCom <- dplyr::select(radCom,Species,Ranks, ys=Commun, Community)}
    if(method == "Log2") { radCom <- dplyr::select(radCom, Species,Ranks, ys=log2pi, Community)}
    if(method == "Log10") { radCom <- dplyr::select(radCom, Species,Ranks, ys=log10pi, Community)}
    if(method == "Rel.Abun") { radCom <- dplyr::select(radCom, Species,Ranks, ys=Re.Ab, Community)}
    
      
    if(name == FALSE){
      if(mixing == FALSE){
        pr <- ggplot2::ggplot(radCom, ggplot2::aes(x = Ranks, y = ys, color = Community)) + geom_point() + facet_wrap(~Community, ncol=2) +
                       geom_line()  + theme_bw() + xlab('Rank Abundance') + ylab(paste0(method)) + 
                       theme(legend.position = "none", axis.text=element_text(size=7), axis.title=element_text(size=7))+
                       scale_fill_brewer(palette="Set1")
      } else{ 
                       
                    pr <- ggplot2::ggplot(radCom, ggplot2::aes(x = Ranks, y = ys, fill = Community, color = Community)) + geom_point() +
                         geom_line()  + theme_bw() + xlab('Rank Abundance') + ylab(paste0(method)) + 
                         theme(legend.position = "right", axis.text=element_text(size=7),
                               axis.title=element_text(size=7)) +
                      scale_color_brewer(palette="Set1")
                       
                       }
    } else{
      
      if(mixing == FALSE){ 
        pr <- ggplot2::ggplot(radCom, ggplot2::aes(x = Ranks, y = ys)) + geom_point() + facet_wrap(vars(Community), ncol=2) +
          geom_line()  + theme_bw() + xlab('Rank Abundance') + ylab(paste0(method)) + 
          theme(legend.position = "none", axis.text=element_text(size=7),
                axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
          geom_text(aes(label=Species),hjust=-0.5, vjust=-0.5) 
      }else{ 
        pr <- ggplot2::ggplot(radCom, ggplot2::aes(x = Ranks, y = ys, fill = Community, color=Community)) + geom_point() + 
          geom_line()  + theme_bw() + xlab('Rank Abundance') + ylab(paste0(method)) + 
          theme(legend.position = "right", axis.text=element_text(size=7),
                axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
            geom_text(aes(label=Species),hjust=-0.5, vjust=-0.5) +
          scale_color_brewer(palette="Set1")
          }
    }
  }
  print(pr)
}


  
 rd <- rad.fit(x, name = T, plot=T, method = 'Rel.Abun', mixing = T)

