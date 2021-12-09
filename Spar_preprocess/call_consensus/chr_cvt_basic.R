## convert yeast chr from numbers to roman numbers 

## the first column have the chr 
chr_num2rom <- function(tb)
{
	chrs <- c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
	for (i in seq(1,16))
	{
		tb [tb[,1] == i,1] <- chrs[i]
	}
	return (tb)
}

## convert from latin to numbers
chr_rom2num <- function(tb)
{
	chrs <- c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
	count <- 1
	for (chr in chrs)
	{
		tb[tb[,1] == chr,1] <- count 
		count <- count +1
	}
	return (tb)
}