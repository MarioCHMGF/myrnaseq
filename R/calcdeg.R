#' A function to perform RNA seq DEG analysis with up to four additional covariates and results annotation
#'
#' @param dgeObjx DGEList from previous analysis
#' @param outcome outcome from sampleinfo
#' @param cov1 covariate in sampleinfo default is NULL
#' @param cov2 covariate in sampleinfo  default is NULL
#' @param cov3 covariate in sampleinfo  default is NULL
#' @param cov4 covariate in sampleinfo  default is NULL
#' @return annotated results
#' @import edgeR
#' @import limma
#' @import gplots
#' @import RColorBrewer
#' @import pheatmap
#' @import dendextend
#' @import org.Hs.eg.db
#'
#' @export

calcdeg = function(dgeObjx, outcome, cov1=NULL, cov2=NULL, cov3=NULL, cov4=NULL){
		if(is.null(cov1)){
		design <- model.matrix(~ outcome)
		y <- voom(dgeObjx, design) 
		prefit <- lmFit(y, design)
		fit <- eBayes(prefit)
		results <- as.data.frame(topTable(fit, coef=2, confint=T, n = Inf))
		ann <- AnnotationDbi::select(org.Hs.eg.db,keys=rownames(results),columns=c("ENTREZID","SYMBOL","GENENAME"))
		results.annotated <- cbind(ann, results)
		write.csv(results.annotated, file="my_unadjusted_results.csv", row.names=FALSE)
		return(results.annotated)
		}
		else
		{
		if(!is.null(cov1)){
		design <- model.matrix(~ outcome + cov1)
		y <- voom(dgeObjx, design) 
		prefit <- lmFit(y, design)
		fit <- eBayes(prefit)
		results <- as.data.frame(topTable(fit, coef=2, confint=T, n = Inf))
		ann <- AnnotationDbi::select(org.Hs.eg.db,keys=rownames(results),columns=c("ENTREZID","SYMBOL","GENENAME"))
		results.annotated <- cbind(ann, results)
		write.csv(results.annotated, file="my_adjusted_results.csv", row.names=FALSE)
		return(results.annotated)
		}
		if(!is.null(c(cov1, cov2))){
			design <- model.matrix(~ outcome + cov1 + cov2)
			y <- voom(dgeObjx, design) 
			prefit <- lmFit(y, design)
			fit <- eBayes(prefit)
			results <- as.data.frame(topTable(fit, coef=2, confint=T, n = Inf))
			ann <- AnnotationDbi::select(org.Hs.eg.db,keys=rownames(results),columns=c("ENTREZID","SYMBOL","GENENAME"))
			results.annotated <- cbind(ann, results)
			write.csv(results.annotated, file="my_adjusted_results.csv", row.names=FALSE)
			return(results.annotated)
			}
				if(!is.null(c(cov1, cov2, cov3))){
				design <- model.matrix(~ outcome + cov1 + cov2 + cov3)
				y <- voom(dgeObjx, design) 
				prefit <- lmFit(y, design)
				fit <- eBayes(prefit)
				results <- as.data.frame(topTable(fit, coef=2, confint=T, n = Inf))
				ann <- AnnotationDbi::select(org.Hs.eg.db,keys=rownames(results),columns=c("ENTREZID","SYMBOL","GENENAME"))
				results.annotated <- cbind(ann, results)
				write.csv(results.annotated, file="my_adjusted_results.csv", row.names=FALSE)
				return(results.annotated)
				}
					if(!is.null(c(cov1, cov2, cov3, cov4))){
					design <- model.matrix(~ outcome + cov1 + cov2 + cov3 + cov4)
					y <- voom(dgeObjx, design) 
					prefit <- lmFit(y, design)
					fit <- eBayes(prefit)
					results <- as.data.frame(topTable(fit, coef=2, confint=T, n = Inf))
					ann <- AnnotationDbi::select(org.Hs.eg.db,keys=rownames(results),columns=c("ENTREZID","SYMBOL","GENENAME"))
					results.annotated <- cbind(ann, results)
					write.csv(results.annotated, file="my_adjusted_results.csv", row.names=FALSE)
					return(results.annotated)
					}
				}
			}
