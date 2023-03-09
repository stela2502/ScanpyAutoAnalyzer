```R
library(SingleR)
library(scuttle)
library(cellexalvrR)
library(SingleCellExperiment)
library(hdf5r)

data = as_cellexalvrR( 'input.h5ad', specie='mouse') # specie is not important
data = SingleCellExperiment( data@data )

toSparse <- function(file){
	message("reading expression data")
	x= file[['X']][['data']][]
	i= file[['X']][['indices']][]
	j= rep(0, length(x))
	indptr = file[['X']][['indptr']][]
	last = 1
	for ( a in 2: length(indptr) ) {
		j[(indptr[a-1]+1):(indptr[a]+1)] = last
		last = last+1
	}
	j = j [-length(j)]
	m = Matrix::sparseMatrix( i = i+1, j=j, x=x)
	meta.data = H5Anno2df( file, 'obs')
	annotation = H5Anno2df( file,'var')
	rownames(m) = annotation[,'_index']
	colnames(m) = meta.data[,'_index']
	m
}

setMethod('H5Anno2df', signature = c ('H5File'),
	definition = function (x, slotName, namecol=NULL, onlyStrings=FALSE ) {
		OK = NULL;
		for (n in names( x[[slotName]])) {
			if ( is(x[[paste(sep="/",slotName,n)]], 'H5Group') ){
				OK= c( OK, FALSE )
			}else {
				OK = c(OK, TRUE)
			}
		}
		obs = data.frame(lapply(names(x[[slotName]])[OK], function(n) { x[[paste(sep="/",slotName,n)]][] } ))
		colnames( obs ) = names(x[[slotName]])[OK]

		obs2 = data.frame()
		try( {
		obs2 = data.frame(lapply(names(x[[slotName]])[!OK], function(n) {
			if ( n == "__categories"){
				return (NULL)
			}
			factor ( 
				x[[paste(sep="/",slotName,n)]][['codes']][], 
				labels= x[[paste(sep="/",slotName,n)]][['categories']][] 
				)
			 } 
		))
		})
		if ( nrow( obs2) == nrow(obs) ){
			colnames( obs2) = names(x[[slotName]])[!OK]
			obs = cbind ( obs, obs2)
		}

		col_uniq= NULL
		for( n in colnames(obs) ) {
			if ( all(obs[,n] =="") ){
				obs[,n] = NULL
			}else {
				col_uniq = c( col_uniq, length(unique(obs[,n]))) 
			}
		}
		names(col_uniq) = colnames( obs )
		## most likely cell names column
		if ( ! is.na(match('_index', colnames( obs )))){
			namecol = '_index'
		}
		if ( ! is.null(namecol )) {
			rownames(obs) =  as.vector(obs[, namecol]) 
		}else {
			## now I need to check for strings...

			OK = unlist(lapply( colnames(obs) , function(id) {
				suppressWarnings({a= which( is.na(as.numeric(as.vector(obs[,id])))==T) })
				## Strings only
				if ( length(a) > 0) {
					length(unique(as.vector(obs[a, id])))
				}else {
					0
				}
			}))

			names(OK) = colnames(obs)
			rownames(obs) =  make.names ( #function definition in file 'as_cellexalvrR.R'
					as.vector(obs[, names(OK)[which(OK == max(OK))[1]]]) )
		}
		## now lets double check the factors
		## x[[slotName]][['__categories']][[id]][]
		try({
			if ( ! is.na( match ('__categories', names(x[[slotName]] )) ) ){

			for (catName in names(x[[slotName]][['__categories']]) ) {
				## if we have strings in the category info
				levelNames = x[[slotName]][['__categories']][[catName]][]
				if ( ! all(!  is.na(suppressWarnings({as.numeric( levelNames )})))){
					## replace the data with an R Factor
					if ( min(obs[,catName]) == -1 ){ ## them would be the NA in R
						# Just add one more category - not perfect, but better than breaking
						levelNames = c( 'NA', levelNames)
						obs[,catName] = factor( levelNames[obs[,catName]+2] )
					}
					else {
						obs[,catName] = factor( levelNames[obs[,catName]+1] )
					}
				}
			}
			}
		})
		if ( onlyStrings ) {
			for( i in 1:length(col_uniq) ) {
				if ( col_uniq[i] == 0 ){ # all entries convertable to numeric
					obs[,i] = NULL
				}
			}
		}

		if ( is.na(match("_index", colnames(obs))) ){
			obs$'_index' = rownames(obs)
		}

		obs
	} 
)