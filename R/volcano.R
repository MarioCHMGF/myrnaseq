#' A function to draw volcano plots
#'
#' @param myresults annotated results from calcdeg
#' @param Q if TRUE take Q values, default is FALSE
#' @param top nu,ber of top genes to annotate, default is 5
#' @param xl lower x axis value, must be negative, default is NULL
#' @param xu upper x axis value, must be positive, default is NULL
#' @return volcanoplot
#' @import ggplot2
#' @import ggrepel
#'
#' @export

volcano = function(myresults, Q=FALSE, top=5, xl=NULL, xu=NULL){
	if(!is.null(xl) && Q==TRUE){
		deg <- myresults
		deg$SIG <- NA
		deg$SIG[deg$adj.P.Val > 0.05] <- "Q > 0.05"
		deg$SIG[deg$adj.P.Val <= 0.05] <- "Q < 0.05"
		deg$LOG <- -log10(deg$adj.P.Val)

		deg$LABEL <- NA
		Top_Hits = head(arrange(deg, adj.P.Val),top)
		deg$LABEL <- if_else(deg$SYMBOL %in% Top_Hits$SYMBOL, deg$SYMBOL, NA)

		plotv <- ggplot(deg, aes(x=logFC, y=LOG, color=SIG, label=LABEL)) +
		geom_point(size=3, alpha = 3/10)+
		geom_vline(xintercept=c(2,-2), color="red", size=0.2,linetype="dashed") +
		geom_vline(xintercept=c(1,-1), color="red", size=0.2,linetype="dashed") +
		geom_hline(yintercept=-log10(0.05), color="red", size=0.2,linetype="dashed") +
		scale_color_manual(values=c('tomato1','grey45')) + xlim(xl, xu) +
		labs(title="Expression",x="Loget", y=expression(-Log[10]~(q-value))) +
		geom_text_repel(color = "black", box.padding = 3, max.overlaps = Inf)
		return(plotv)
		}
	if(is.null(xl) && Q==TRUE){
		deg <- myresults
		deg$SIG <- NA
		deg$SIG[deg$adj.P.Val > 0.05] <- "Q > 0.05"
		deg$SIG[deg$adj.P.Val <= 0.05] <- "Q < 0.05"
		deg$LOG <- -log10(deg$adj.P.Val)

		deg$LABEL <- NA
		Top_Hits = head(arrange(deg, adj.P.Val),top)
		deg$LABEL <- if_else(deg$SYMBOL %in% Top_Hits$SYMBOL, deg$SYMBOL, NA)

		plotv <- ggplot(deg, aes(x=logFC, y=LOG, color=SIG, label=LABEL)) +
		geom_point(size=3, alpha = 3/10)+
		geom_vline(xintercept=c(2,-2), color="red", size=0.2,linetype="dashed") +
		geom_vline(xintercept=c(1,-1), color="red", size=0.2,linetype="dashed") +
		geom_hline(yintercept=-log10(0.05), color="red", size=0.2,linetype="dashed") +
		scale_color_manual(values=c('tomato1','grey45')) +
		labs(title="Expression",x="Loget", y=expression(-Log[10]~(q-value))) +
		geom_text_repel(color = "black", box.padding = 3, max.overlaps = Inf)
		return(plotv)
		}
	if(!is.null(xl) && Q==FALSE){
		deg <- myresults
		deg$SIG <- NA
		deg$SIG[deg$P.Value > 0.05] <- "P > 0.05"
		deg$SIG[deg$P.Value <= 0.05] <- "P < 0.05"
		deg$LOG <- -log10(deg$P.Value)

		deg$LABEL <- NA
		Top_Hits = head(arrange(deg, P.Value),top)
		deg$LABEL <- if_else(deg$SYMBOL %in% Top_Hits$SYMBOL, deg$SYMBOL, NA)

		plotv <- ggplot(deg, aes(x=logFC, y=LOG, color=SIG, label=LABEL)) +
		geom_point(size=3, alpha = 3/10)+
		geom_vline(xintercept=c(2,-2), color="red", size=0.2,linetype="dashed") +
		geom_vline(xintercept=c(1,-1), color="red", size=0.2,linetype="dashed") +
		geom_hline(yintercept=-log10(0.05), color="red", size=0.2,linetype="dashed") +
		scale_color_manual(values=c('tomato1','grey45')) + xlim(xl, xu) +
		labs(title="Expression",x="Loget", y=expression(-Log[10]~(p-value))) +
		geom_text_repel(color = "black", box.padding = 3, max.overlaps = Inf)
		return(plotv)
		}
	if(is.null(xl) && Q==FALSE){
		deg <- myresults
		deg$SIG <- NA
		deg$SIG[deg$P.Value > 0.05] <- "P > 0.05"
		deg$SIG[deg$P.Value <= 0.05] <- "P < 0.05"
		deg$LOG <- -log10(deg$P.Value)

		deg$LABEL <- NA
		Top_Hits = head(arrange(deg, P.Value),top)
		deg$LABEL <- if_else(deg$SYMBOL %in% Top_Hits$SYMBOL, deg$SYMBOL, NA)

		plotv <- ggplot(deg, aes(x=logFC, y=LOG, color=SIG, label=LABEL)) +
		geom_point(size=3, alpha = 3/10)+
		geom_vline(xintercept=c(2,-2), color="red", size=0.2,linetype="dashed") +
		geom_vline(xintercept=c(1,-1), color="red", size=0.2,linetype="dashed") +
		geom_hline(yintercept=-log10(0.05), color="red", size=0.2,linetype="dashed") +
		scale_color_manual(values=c('tomato1','grey45')) +
		labs(title="Expression",x="Loget", y=expression(-Log[10]~(p-value))) +
		geom_text_repel(color = "black", box.padding = 3, max.overlaps = Inf)
		return(plotv)
		}
	}
