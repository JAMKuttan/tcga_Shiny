library(ggplot2)
library(shiny)
library(shinyjs)
library("survival")
library("survminer")

get.filt<-function(var){
  if(var$operation=="Pancan"){
    return(
      list(
        operation=var$operation,
        gene=var$gene,
        cancer.list=var$cancerType
      )
    )
  }else{
    return(
      list(
        operation=var$operation,
        gene=var$gene2,
        cancer=var$cancerType2,
        cutoff=var$cut.off
      )
    )
  }
}

shinyServer(
  function(input,output,session){
    
    vars<-eventReactive(input$go,{
      get.filt(input)
    })
    
    output$myPlot <- renderPlot({
      if(vars()$operation=="Pancan"){

        gene<-as.character(vars()$gene)
        cancer.list<-as.character(vars()$cancer.list)
        
        data <- data.frame(value=as.numeric(),type=as.character(),cancer=as.character(),pval=as.numeric())
        pval <- NULL
        for(i in cancer.list){
          if(!file.exists(paste0("../Expression/out_genes/",gene,"/",i,".txt"))){
            cat("File not present\n")
            system(paste0("mkdir ../Expression/out_genes/",gene))
            cat("Made directory/n")
            system(paste("ls ../Expression/*txt | while read line ; do fname=$(basename $line); head -1 $line > ../Expression/out_genes/",
                         gene,"/","$fname; grep -wi ",gene," $line >> ../Expression/out_genes/",gene,"/","$fname;done",sep=""))
            cat("creating data\n")
            
          }
          temp <- as.data.frame(t(read.table(paste("../Expression/out_genes/",gene,"/",i,".txt",sep=""),header=T,sep="\t",row.names=1)))
          temp$type <- as.factor(ifelse(grepl("11",row.names(temp)),"Normal","Tumor"))
          temp$cancer <- rep(i,nrow(temp))
          if(ncol(temp)<3 | length(unique(temp$type))<2){ 
            next 
          }else{
            colnames(temp) <- c("value","type","cancer")
            print(i)
            pval <- c(pval,signif(kruskal.test(value~type,temp)$p.value,2))
            data <- rbind(data,temp)
          }
        }
        
        if(nrow(data)==0){
          cat("Not enough data for this cancer type\n")
        }else{
          data$value <- 2^data$value
          p <- ggplot(data = data, aes(x=type, y=value)) + geom_boxplot(aes(fill=type),outlier.size = 0.2,outlier.alpha=0.8) + scale_y_continuous(limits = c(0, 2000)) + geom_text(data=data.frame(x=1.5,y=1800,labels=pval,cancer=cancer.list),aes(x,y,label=labels),size=3 , inherit.aes=FALSE)
          #cat("printing",gene,"\n")
          print(p + facet_wrap( ~ cancer, scales="free") + ggtitle(gene) + theme(plot.title = element_text(hjust = 0.5)))

        }
      }else if(vars()$operation=="Survival"){
        #Define Cancer & Gene
        cat("Survival plot/n")
        cancer <- vars()$cancer
        gene <- vars()$gene
        cat(paste0(cancer," ",gene,"\n"))
        
        #Import Expression and Clinical Data
        cat(paste0("going to read expression ",cancer))
        exp <- read.table(paste0("../Expression/",cancer,".txt"),header=T,row.names=1,check.names=F) #Instead of importing entire expression matrix should import single gene.
        # exp <- read.table(paste("Expression/",cancer,".txt",sep=""),header=T,row.names=1,check.names=F)
        cat("going to read clinical data")
        clinical <- read.table(paste("../Clinical/",cancer,".txt",sep=""),header=T,row.names=1,sep="\t")
        
        
        #Process Clinical
        # get the columns that contain data we can use: days to death, new tumor event, last day contact to....
        ind_keep <- grep("days_to_new_tumor_event_after_initial_treatment",colnames(clinical))
        
        # this is a bit tedious, since there are numerous follow ups, let's collapse them together and keep the first value (the higher one) if more than one is available
        new_tum <- as.matrix(clinical[,ind_keep])
        new_tum_collapsed <- c()
        for (i in 1:dim(new_tum)[1]){
          if(sum(is.na(new_tum[i,])) < dim(new_tum)[2]){
            m <- min(new_tum[i,],na.rm=T)
            new_tum_collapsed <- c(new_tum_collapsed,m)
          } else {
            new_tum_collapsed <- c(new_tum_collapsed,"NA")
          }
        }
        
        # do the same to death
        ind_keep <- grep("days_to_death",colnames(clinical))
        death <- as.matrix(clinical[,ind_keep])
        death_collapsed <- c()
        for (i in 1:dim(death)[1]){
          if(sum(is.na(death[i,])) < dim(death)[2]){
            m <- max(death[i,],na.rm=T)
            death_collapsed <- c(death_collapsed,m)
          } else {
            death_collapsed <- c(death_collapsed,"NA")
          }
        }
        
        # and days last follow up here we take the most recent which is the max number
        ind_keep <- grep("days_to_last_followup",colnames(clinical))
        fl <- as.matrix(clinical[,ind_keep])
        fl_collapsed <- c()
        for (i in 1:dim(fl)[1]){
          if(sum(is.na(fl[i,])) < dim(fl)[2]){
            m <- max(fl[i,],na.rm=T)
            fl_collapsed <- c(fl_collapsed,m)
          } else {
            fl_collapsed <- c(fl_collapsed,"NA")
          }
        }
        
        # and put everything together
        all_clin <- data.frame(new_tum_collapsed,death_collapsed,fl_collapsed)
        colnames(all_clin) <- c("new_tumor_days", "death_days", "followUp_days")
        
        # create vector with time to new tumor containing data to censor for new_tumor
        all_clin$new_time <- c()
        for (i in 1:length(as.numeric(as.character(all_clin$new_tumor_days)))){
          all_clin$new_time[i] <- ifelse(is.na(as.numeric(as.character(all_clin$new_tumor_days))[i]),
                                         as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$new_tumor_days))[i])
        }
        
        # create vector time to death containing values to censor for death
        all_clin$new_death <- c()
        for (i in 1:length(as.numeric(as.character(all_clin$death_days)))){
          all_clin$new_death[i] <- ifelse(is.na(as.numeric(as.character(all_clin$death_days))[i]),
                                          as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$death_days))[i])
        }
        
        all_clin$death_event <- ifelse(clinical$vital_status == "LIVING", 0,1)
        row.names(all_clin) <- row.names(clinical)
        
        cutoff <- as.numeric(vars()$cutoff)
        #Setup Data
        ### Deprecated Code ###
        # data <- all_clin[names[names %in% row.names(all_clin)],]#[which(substr(colnames(exp),14,14) == "0"),]
        #all_clin$event_rna <- apply(exp[gene,], 1, function(x) ifelse(abs(x) > median(as.numeric(exp[gene,])),1,0))
        # gene.exp <- apply(gene.exp,2,function(x) ifelse(x > quantile(as.numeric(x),1-cutoff),1,0))
        # gene.exp <- t(as.data.frame(sapply(by(gene.exp,rownames(gene.exp),colSums),identity)))
        colnames(exp) <- gsub("01A","01",colnames(exp))
        gene.exp <- t(exp[gene,])
        colnames(gene.exp) <- "gene"
        data <- merge(all_clin,gene.exp,by="row.names")
        data$event_rna <- ifelse(data$gene > quantile(as.numeric(data$gene),1-cutoff),1,0)
        data <- data[which(data$gene > quantile(as.numeric(data$gene),1-cutoff) | data$gene < quantile(as.numeric(data$gene),cutoff)),]
        # run survival analysis
        s <- surv_fit(Surv(as.numeric(as.character(data$new_death)),data$death_event)~event_rna,data)
        p<-ggsurvplot(s,conf.int=F,xlab="Time(days)",pval=T,pval.coord=c(1000,0.3),risk.table=TRUE,legend.title="",legend=c(0.8,0.8),legend.labs = c( paste(gene,"LOW",sep=" "), paste(gene,"HIGH",sep=" ")))
        print(p+ggtitle(paste(cancer, gene,"cutoff: ", cutoff)))
      }
      
    })
    
    
  }
)
