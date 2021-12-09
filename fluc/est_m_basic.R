## basic functons for estimate m

library(ggplot2)
library("rsalvador")
library(dplyr)

## this function will apply to one fluctuation experiment (multiple strains, one experiment)
est_m <- function(mut_tb, strains)
{
    mut_m_vec <- NULL
    mut_m_ci_vec <- list()
  
    for(strain in strains)
    {
      strain_mut <- mut_tb [mut_tb[,1] == strain,4]
      strain_mut  <- strain_mut[!is.na(strain_mut)]
      #cells <- rep(cell_count, length(strain_mut))
    
      mut_m <- newton.LD(strain_mut)
      m_ci <- list(confint.LD(strain_mut))
      mut_m_vec <- c(mut_m_vec, mut_m )
    
      mut_m_ci_vec <- append(mut_m_ci_vec , m_ci )
    }
    
    return (list(mut_m_vec, mut_m_ci_vec))
}

## this function applies to the data which has multiple experiments combined (multiple strains, multiple experiments) 
## will use column 2 "Expt" to distinguish different exp
est_m_multi <- function (mut_file)
{
  mut_tb <- read.table(mut_file,sep=",", header=TRUE, stringsAsFactors = FALSE)
  
  ## estimate m for each experiment (ignore confidence interval for now)
  mut_tb %>% group_by(Strain, Expt) %>% summarise(newton.LD(mut_colonies)) %>% data.frame -> exp_m 
  
  names(exp_m) <- c("Strain", "Expt","m")
  
  return (exp_m)
}

## calculate mutation rate from mut table and cell count table
calc_mut_rate_all <- function (mut_file, count_file)
{
  m_tb <- est_m_multi(mut_file)
  
  cell_count <- read.table(count_file, sep=",",header= TRUE, stringsAsFactors = FALSE)
  
  ## merge m table and cell count table 
  
  tb <- merge(cell_count, m_tb, by=c("Strain", "Expt"))
  
  tb[,"m_rate"] <- tb[,"m"]/tb[,"cell_count"]
  
  return (tb)
}

## tb is the output from calc_mut_rate_all
## add color_col (the column for color)
plot_mut_rate <- function(tb, color_col=NULL, shape_col=NULL)
{
  
  ## first get the mean for m_rate for each strain 
  tb %>% group_by(Strain) %>% summarise(mean(m_rate) )  %>%  data.frame -> m_rate
  
  ## order strain by the mut rate 
  strain_order <- as.character (m_rate[order(m_rate[,2]),1])
  
  ## should create a new factor obj and assign levels, 
  ## should not just change levels of currect factor (when plotting, it only changes x labels without changing y with it)
  tb[,"Strain"] <- factor(tb[,"Strain"], levels= strain_order)
  
  if(is.null(color_col))
  {
    p <- ggplot (tb, aes(Strain,m_rate)) + geom_boxplot() + geom_point() + theme(axis.text=element_text(size=9))
  }
  else
  {
    if(is.null(shape_col))
    {
      p <- ggplot (tb, aes_string("Strain","m_rate")) + geom_boxplot(aes_string(color= color_col)) +
        geom_point(aes_string(color= color_col), size=1.8) + theme(axis.text=element_text(size=9))
    }
    else
    {
      p <- ggplot (tb, aes_string("Strain","m_rate")) + geom_boxplot(aes_string(color= color_col)) +
        geom_point(aes_string(color= color_col, shape= shape_col), size=2) +
        theme(axis.text=element_text(size=9))
    }
    
  }
  
  return (p)
}


## did not order by mutation rate, but by clade 
plot_mut_rate2 <- function(tb, color_col=NULL, shape_col=NULL)
{
  
  ## first get the mean for m_rate for each strain 
  tb %>% group_by(Strain) %>% summarise(mean(m_rate) )  %>%  data.frame -> m_rate
  names(m_rate) <- c("Strain", "strain_m_rate")
  
  ## also calculate mean mut rate for a clade  
  tb %>% group_by(Clades) %>% summarise(mean(m_rate) )  %>%  data.frame -> m_rate_clade
  names(m_rate_clade) <- c("Clades","clades_m_rate")
   
  tb_subcol <- tb [,c("Strain","Clades")]
  
  mean_mrate_tb <- merge(m_rate, tb_subcol, by="Strain")
  mean_mrate_tb2 <- merge(m_rate_clade, mean_mrate_tb, by="Clades")
  
  strain_order <- unique(mean_mrate_tb2 [order(mean_mrate_tb2[,"clades_m_rate"], mean_mrate_tb2[,"strain_m_rate"]),"Strain"])
  
  
  ## should create a new factor obj and assign levels, 
  ## should not just change levels of currect factor (when plotting, it only changes x labels without changing y with it)
  tb[,"Strain"] <- factor(tb[,"Strain"], levels= strain_order)
  
  if(is.null(color_col))
  {
    p <- ggplot (tb, aes(Strain,m_rate)) + geom_boxplot() + geom_point() + theme(axis.text=element_text(size=9))
  }
  else
  {
    if(is.null(shape_col))
    {
      p <- ggplot (tb, aes_string("Strain","m_rate")) + geom_boxplot(aes_string(color= color_col)) +
        geom_point(aes_string(color= color_col), size=1.8) + theme(axis.text=element_text(size=9))
    }
    else
    {
      p <- ggplot (tb, aes_string(x="Strain",y="m_rate", fill="Clades")) + geom_boxplot(aes_string(color= color_col)) +
        geom_point(aes_string(color= color_col, shape= shape_col), size=2) +
        theme(axis.text=element_text(size=9))
    }
    
  }
  
  return (p)
}

