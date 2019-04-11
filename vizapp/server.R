library(ggplot2)
library(shiny)
library(shinyjs)
library("survival")
library("survminer")

cat("reading in gene list...")
all.genes<-read.table("tcga_genes.txt",sep = "\t",header = FALSE)
all.genes<-as.character(all.genes[,1])
cat("done\n")
get.filt<-function(var){
  if(var$operation=="Boxplot"){
    return(
      list(
        operation=var$operation,
        gene=var$gene,
        cancer.list=var$cancerType
      )
    )
  }else{
    if(var$operation2 == "Manual"){
      return(
        list(
          operation=var$operation,
          operation2=var$operation2,
          gene=var$gene2,
          cancer=var$cancerType2,
          cutoff=var$cut.off
        )
      )
    }else if(var$operation2 == "Automatic"){
      return(
        list(
          operation=var$operation,
          operation2=var$operation2,
          gene=var$gene2,
          cancer=var$cancerType2
        )
      )
    }
  }
}

shinyServer(
  function(input,output,session){
    
    vars<-eventReactive(input$go,{
      get.filt(input)
    })
    
    output$myPlot <- renderPlot({
      #PANCAN BOX PLOT
      withProgress(message = 'Making plot', value = 0,{
        if(vars()$operation=="Boxplot"){
          
          gene<-toupper(gsub("\\s", "", as.character(vars()$gene)))
          cancer.list<-as.character(vars()$cancer.list)
          if(!(gene %in% all.genes)){
            showNotification(paste0("SORRY! ", gene, " is not a valid gene!"), duration = 8, type = "error")
          }else{
            data <- data.frame(value=as.numeric(),type=as.character(),cancer=as.character(),pval=as.numeric())
            pval <- NULL
            count <- 1
            n = length(cancer.list) + 1
            for(i in cancer.list){
              incProgress(1/n, detail = paste0("collating data ", count, " of ", length(cancer.list)))
              if(!file.exists(paste0("../output/Expression/out_genes/",gene,"/",i,".txt"))){
                cat("File not present\n")
                if(!dir.exists(paste0("../output/Expression/out_genes/",gene))){
                  system(paste0("mkdir ../output/Expression/out_genes/",gene))
                  cat("Made directory\n")
                }
                system(paste0("ls ../output/Expression/*",i,"*txt | while read line ; do fname=$(basename $line); head -1 $line > ../output/Expression/out_genes/",
                              gene,"/","$fname; grep -wi ",gene," $line >> ../output/Expression/out_genes/",gene,"/","$fname;done"))
                cat("creating data\n")
                
              }
              temp <- as.data.frame(t(read.table(paste("../output/Expression/out_genes/",gene,"/",i,".txt",sep=""),header=T,sep="\t",row.names=1)))
              temp$type <- as.factor(ifelse(grepl("11",row.names(temp)),"Normal","Tumor"))
              temp$cancer <- rep(i,nrow(temp))
              if(ncol(temp)<3 | length(unique(temp$type))<2){
                if(length(cancer.list) > 1){
                  showNotification(paste0("WARNING! Project ", i, " does not have normal tissue data! Skipping..."), duration = 8, type = "warning")
                  next 
                }else{
                  showNotification(paste0("SORRY! Project ", i, " does not have normal tissue data! Choose another project."), duration = 8, type = "error")
                  break
                }
                
              }else{
                colnames(temp) <- c("value","type","cancer")
                print(i)
                pval <- c(pval,signif(kruskal.test(value~type,temp)$p.value,2))
                data <- rbind(data,temp)
              }
              count <- count + 1
            }
            
            if(nrow(data)==0){
              return()
            }else{
              data$value <- 2^data$value
              p <- ggplot(data = data, aes(x=type, y=value)) + geom_boxplot(aes(fill=type),outlier.size = 0.2,outlier.alpha=0.8) + scale_y_continuous(limits = c(0, 2000)) + 
                geom_text(data=data.frame(x=1.5,y=1800,labels=pval,cancer=unique(data$cancer)),aes(x,y,label=labels),size=3 , inherit.aes=FALSE) +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black"))
              #cat("printing",gene,"\n")
              print(p + facet_wrap( ~ cancer, scales="free") + ggtitle(gene) + theme(plot.title = element_text(hjust = 0.5)))
              incProgress(1/n, detail = paste0("printing plot!"))
            }
          }
          #SURVIVAL PLOTS  
        }else if(vars()$operation=="Survival"){
          #Define Cancer & Gene
          oldw <- getOption("warn")
          options(warn = -1)
          cat("Survival plot\n")
          cancer <- vars()$cancer
          gene<-toupper(gsub("\\s", "", as.character(vars()$gene)))
          cat(paste0(cancer," ",gene,"\n"))
          if(!(gene %in% all.genes)){
            showNotification(paste0("SORRY! ", gene, " is not a valid gene!"), duration = 8, type = "error")
          }else{
            n <- 4
            incProgress(1/n, detail = "Getting expression data")
            #Import Expression and Clinical Data
            cat(paste0("going to read expression ",cancer,"\n"))
            if(!file.exists(paste0("../output/Expression/out_genes/",gene,"/",cancer,".txt"))){
              cat("File not present\n")
              if(!dir.exists(paste0("../output/Expression/out_genes/",gene))){
                system(paste0("mkdir ../output/Expression/out_genes/",gene))
                cat("Made directory\n")
              }
              system(paste0("ls ../output/Expression/*",cancer,"*txt | while read line ; do fname=$(basename $line); head -1 $line > ../output/Expression/out_genes/",
                            gene,"/","$fname; grep -wi ",gene," $line >> ../output/Expression/out_genes/",gene,"/","$fname;done"))
              cat("creating data\n")
              
            }
            exp <- read.table(paste0("../output/Expression/out_genes/",gene,"/",cancer,".txt"),header=T,row.names=1,check.names=F) 
            #exp <- read.table(paste0("../output/Expression/",cancer,".txt"),header=T,row.names=1,check.names=F) #Instead of importing entire expression matrix should import single gene.
            cat("going to read clinical data\n")
            incProgress(1/n, detail = "Getting clinical data")
            clinical <- read.table(paste("../output/Clinical/",cancer,".txt",sep=""),header=T,row.names=1,sep="\t")
            
            
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
            colnames(exp) <- gsub("01A","01",colnames(exp))
            gene.exp <- t(exp[gene,])
            colnames(gene.exp) <- "gene"
            
            #Fork for Manual vs Automatic selection
            if(vars()$operation2 == "Manual"){
              cutoff <- as.numeric(vars()$cutoff)
              #Setup Data

              colnames(exp) <- gsub("01A","01",colnames(exp))
              gene.exp <- t(exp[gene,])
              colnames(gene.exp) <- "gene"
              data <- merge(all_clin,gene.exp,by="row.names")
              data$event_rna <- ifelse(data$gene > quantile(as.numeric(data$gene),1-cutoff),1,0)
              data <- data[which(data$gene > quantile(as.numeric(data$gene),1-cutoff) | data$gene < quantile(as.numeric(data$gene),cutoff)),]
              # run survival analysis
              s <- surv_fit(Surv(as.numeric(as.character(data$new_death)),data$death_event)~event_rna,data)
              incProgress(1/n, detail = "Printing plot")
              p<-ggsurvplot(s,conf.int=F,xlab="Time(days)",pval=T,pval.coord=c(1000,0.3),risk.table=TRUE,legend.title="",legend=c(0.8,0.8),legend.labs = c( paste(gene,"LOW",sep=" "), paste(gene,"HIGH",sep=" ")))
              print(p+ggtitle(paste(cancer, gene,"cutoff: ", cutoff)))
              options(warn = oldw)
            }else{
              cutoff.vals=c(.10, .15, .20, .25, .30, .35, .40, .45, .50)
              incProgress(1/n, detail = "Checking Cutoffs")
              for(cutoff in cutoff.vals){
                cat("\tchecking cutoff value:",cutoff,"\n")  ###add status update here
                data <- merge(all_clin,gene.exp,by="row.names")
                if(is.na(data$gene)[1]){
                  showNotification(paste0("SORRY! No data for gene ", gene,"!"), duration = 8, type = "error")
                  break
                }
                data$event_rna <- ifelse(data$gene > quantile(as.numeric(data$gene),1-cutoff),1,0)
                data <- data[which(data$gene > quantile(as.numeric(data$gene),1-cutoff) | data$gene < quantile(as.numeric(data$gene),cutoff)),]
                # run survival analysis
                if(length(unique(data$event_rna)) < 2){
                  next
                }
                
                cat("done!\n")
                if(cutoff == 0.10){
                  cat("first run\n")
                  best.data <- data
                  best.res <- surv_fit(Surv(as.numeric(as.character(best.data$new_death)),best.data$death_event)~event_rna,best.data)
                  p<-ggsurvplot(best.res,conf.int=F,xlab="Time(days)",pval=TRUE,pval.coord=c(1000,0.3),risk.table=TRUE,legend.title="",legend=c(0.8,0.8),
                                legend.labs = c( paste(gene,"LOW",sep=" "), paste(gene,"HIGH",sep=" ")))
                  best.cutoff <- cutoff
                }else{
                  s <- surv_fit(Surv(as.numeric(as.character(data$new_death)),data$death_event)~event_rna,data)
                  p<-ggsurvplot(s,conf.int=F,xlab="Time(days)",pval=TRUE,pval.coord=c(1000,0.3),risk.table=TRUE,legend.title="",legend=c(0.8,0.8),
                                legend.labs = c( paste(gene,"LOW",sep=" "), paste(gene,"HIGH",sep=" ")))
                  if((as.numeric(surv_pvalue(s)[2]) < as.numeric(surv_pvalue(best.res)[2])) & (p$table$data$n.risk[1] > 10) & (p$table$data$n.risk[6] > 10)){
                    best.data <- data
                    best.res <- surv_fit(Surv(as.numeric(as.character(best.data$new_death)),best.data$death_event)~event_rna,best.data)
                    p<-ggsurvplot(best.res,conf.int=F,xlab="Time(days)",pval=TRUE,pval.coord=c(1000,0.3),risk.table=TRUE,legend.title="",legend=c(0.8,0.8),
                                  legend.labs = c( paste(gene,"LOW",sep=" "), paste(gene,"HIGH",sep=" ")))
                    best.cutoff <- cutoff
                  }
                }  

              }
              incProgress(1/n, detail = "Printing plot")
              p <- ggsurvplot(best.res,conf.int=F,xlab="Time(days)",pval=TRUE,pval.coord=c(1000,0.3),risk.table=TRUE,legend.title="",legend=c(0.8,0.8),
                              legend.labs = c( paste(gene,"LOW",sep=" "), paste(gene,"HIGH",sep=" ")))
              print(p+ggtitle(paste(cancer, "cutoff: ", best.cutoff)))
              
            }
          }
        }
      })
      
    })
    
    plotInput = function(){
      if(vars()$operation=="Boxplot"){
        
        gene<-toupper(gsub("\\s", "", as.character(vars()$gene)))
        cancer.list<-as.character(vars()$cancer.list)
        if(!(gene %in% all.genes)){
          #showNotification(paste0("SORRY! ", gene, " is not a valid gene!"), duration = 8, type = "error")
        }else{
          data <- data.frame(value=as.numeric(),type=as.character(),cancer=as.character(),pval=as.numeric())
          pval <- NULL
          count <- 1
          n = length(cancer.list) + 1
          for(i in cancer.list){
            #incProgress(1/n, detail = paste0("collating data ", count, " of ", length(cancer.list)))
            if(!file.exists(paste0("../output/Expression/out_genes/",gene,"/",i,".txt"))){
              cat("File not present\n")
              if(!dir.exists(paste0("../output/Expression/out_genes/",gene))){
                system(paste0("mkdir ../output/Expression/out_genes/",gene))
                cat("Made directory\n")
              }
              system(paste0("ls ../output/Expression/*",i,"*txt | while read line ; do fname=$(basename $line); head -1 $line > ../output/Expression/out_genes/",
                            gene,"/","$fname; grep -wi ",gene," $line >> ../output/Expression/out_genes/",gene,"/","$fname;done"))
              cat("creating data\n")
              
            }
            temp <- as.data.frame(t(read.table(paste("../output/Expression/out_genes/",gene,"/",i,".txt",sep=""),header=T,sep="\t",row.names=1)))
            temp$type <- as.factor(ifelse(grepl("11",row.names(temp)),"Normal","Tumor"))
            temp$cancer <- rep(i,nrow(temp))
            if(ncol(temp)<3 | length(unique(temp$type))<2){
              if(length(cancer.list) > 1){
                #showNotification(paste0("WARNING! Project ", i, " does not have normal tissue data! Skipping..."), duration = 8, type = "warning")
                next 
              }else{
                #showNotification(paste0("SORRY! Project ", i, " does not have normal tissue data! Choose another project."), duration = 8, type = "error")
                break
              }
              
            }else{
              colnames(temp) <- c("value","type","cancer")
              print(i)
              pval <- c(pval,signif(kruskal.test(value~type,temp)$p.value,2))
              data <- rbind(data,temp)
            }
            count <- count + 1
          }
          
          if(nrow(data)==0){
            return()
          }else{
            data$value <- 2^data$value
            p <- ggplot(data = data, aes(x=type, y=value)) + geom_boxplot(aes(fill=type),outlier.size = 0.2,outlier.alpha=0.8) + scale_y_continuous(limits = c(0, 2000)) + 
              geom_text(data=data.frame(x=1.5,y=1800,labels=pval,cancer=unique(data$cancer)),aes(x,y,label=labels),size=3 , inherit.aes=FALSE) +
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(colour = "black"))
            #cat("printing",gene,"\n")
            print(p + facet_wrap( ~ cancer, scales="free") + ggtitle(gene) + theme(plot.title = element_text(hjust = 0.5)))
            #incProgress(1/n, detail = paste0("printing plot!"))
          }
        }
        #SURVIVAL PLOTS  
      }else if(vars()$operation=="Survival"){
        #Define Cancer & Gene
        oldw <- getOption("warn")
        options(warn = -1)
        cat("Survival plot\n")
        cancer <- vars()$cancer
        gene<-toupper(gsub("\\s", "", as.character(vars()$gene)))
        cat(paste0(cancer," ",gene,"\n"))
        if(!(gene %in% all.genes)){
          #showNotification(paste0("SORRY! ", gene, " is not a valid gene!"), duration = 8, type = "error")
        }else{
          n <- 4
          #incProgress(1/n, detail = "Getting expression data")
          #Import Expression and Clinical Data
          cat(paste0("going to read expression ",cancer,"\n"))
          if(!file.exists(paste0("../output/Expression/out_genes/",gene,"/",cancer,".txt"))){
            cat("File not present\n")
            if(!dir.exists(paste0("../output/Expression/out_genes/",gene))){
              system(paste0("mkdir ../output/Expression/out_genes/",gene))
              cat("Made directory\n")
            }
            system(paste0("ls ../output/Expression/*",cancer,"*txt | while read line ; do fname=$(basename $line); head -1 $line > ../output/Expression/out_genes/",
                          gene,"/","$fname; grep -wi ",gene," $line >> ../output/Expression/out_genes/",gene,"/","$fname;done"))
            cat("creating data\n")
            
          }
          exp <- read.table(paste0("../output/Expression/out_genes/",gene,"/",cancer,".txt"),header=T,row.names=1,check.names=F) 
          #exp <- read.table(paste0("../output/Expression/",cancer,".txt"),header=T,row.names=1,check.names=F) #Instead of importing entire expression matrix should import single gene.
          cat("going to read clinical data\n")
          #incProgress(1/n, detail = "Getting clinical data")
          clinical <- read.table(paste("../output/Clinical/",cancer,".txt",sep=""),header=T,row.names=1,sep="\t")
          
          
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
          colnames(exp) <- gsub("01A","01",colnames(exp))
          gene.exp <- t(exp[gene,])
          colnames(gene.exp) <- "gene"
          
          #Fork for Manual vs Automatic selection
          if(vars()$operation2 == "Manual"){
            cutoff <- as.numeric(vars()$cutoff)
            #Setup Data
            
            colnames(exp) <- gsub("01A","01",colnames(exp))
            gene.exp <- t(exp[gene,])
            colnames(gene.exp) <- "gene"
            data <- merge(all_clin,gene.exp,by="row.names")
            data$event_rna <- ifelse(data$gene > quantile(as.numeric(data$gene),1-cutoff),1,0)
            data <- data[which(data$gene > quantile(as.numeric(data$gene),1-cutoff) | data$gene < quantile(as.numeric(data$gene),cutoff)),]
            # run survival analysis
            s <- surv_fit(Surv(as.numeric(as.character(data$new_death)),data$death_event)~event_rna,data)
            #incProgress(1/n, detail = "Printing plot")
            p<-ggsurvplot(s,conf.int=F,xlab="Time(days)",pval=T,pval.coord=c(1000,0.3),risk.table=TRUE,legend.title="",legend=c(0.8,0.8),legend.labs = c( paste(gene,"LOW",sep=" "), paste(gene,"HIGH",sep=" ")))
            print(p+ggtitle(paste(cancer, gene,"cutoff: ", cutoff)))
            options(warn = oldw)
          }else{
            cutoff.vals=c(.10, .15, .20, .25, .30, .35, .40, .45, .50)
            #incProgress(1/n, detail = "Checking Cutoffs")
            for(cutoff in cutoff.vals){
              cat("\tchecking cutoff value:",cutoff,"\n")  ###add status update here
              data <- merge(all_clin,gene.exp,by="row.names")
              if(is.na(data$gene)[1]){
                showNotification(paste0("SORRY! No data for gene ", gene,"!"), duration = 8, type = "error")
                break
              }
              data$event_rna <- ifelse(data$gene > quantile(as.numeric(data$gene),1-cutoff),1,0)
              data <- data[which(data$gene > quantile(as.numeric(data$gene),1-cutoff) | data$gene < quantile(as.numeric(data$gene),cutoff)),]
              # run survival analysis
              if(length(unique(data$event_rna)) < 2){
                next
              }
              
              cat("done!\n")
              if(cutoff == 0.10){
                cat("first run\n")
                best.data <- data
                best.res <- surv_fit(Surv(as.numeric(as.character(best.data$new_death)),best.data$death_event)~event_rna,best.data)
                p<-ggsurvplot(best.res,conf.int=F,xlab="Time(days)",pval=TRUE,pval.coord=c(1000,0.3),risk.table=TRUE,legend.title="",legend=c(0.8,0.8),
                              legend.labs = c( paste(gene,"LOW",sep=" "), paste(gene,"HIGH",sep=" ")))
                best.cutoff <- cutoff
              }else{
                s <- surv_fit(Surv(as.numeric(as.character(data$new_death)),data$death_event)~event_rna,data)
                p<-ggsurvplot(s,conf.int=F,xlab="Time(days)",pval=TRUE,pval.coord=c(1000,0.3),risk.table=TRUE,legend.title="",legend=c(0.8,0.8),
                              legend.labs = c( paste(gene,"LOW",sep=" "), paste(gene,"HIGH",sep=" ")))
                if((as.numeric(surv_pvalue(s)[2]) < as.numeric(surv_pvalue(best.res)[2])) & (p$table$data$n.risk[1] > 10) & (p$table$data$n.risk[6] > 10)){
                  best.data <- data
                  best.res <- surv_fit(Surv(as.numeric(as.character(best.data$new_death)),best.data$death_event)~event_rna,best.data)
                  p<-ggsurvplot(best.res,conf.int=F,xlab="Time(days)",pval=TRUE,pval.coord=c(1000,0.3),risk.table=TRUE,legend.title="",legend=c(0.8,0.8),
                                legend.labs = c( paste(gene,"LOW",sep=" "), paste(gene,"HIGH",sep=" ")))
                  best.cutoff <- cutoff
                }
              }  
              
            }
            #incProgress(1/n, detail = "Printing plot")
            p <- ggsurvplot(best.res,conf.int=F,xlab="Time(days)",pval=TRUE,pval.coord=c(1000,0.3),risk.table=TRUE,legend.title="",legend=c(0.8,0.8),
                            legend.labs = c( paste(gene,"LOW",sep=" "), paste(gene,"HIGH",sep=" ")))
            print(p+ggtitle(paste(cancer, "cutoff: ", best.cutoff)))
            
          }
        }
      }
    }
    
    output$downloadplot <- downloadHandler(
      filename = function() { paste0("plot", '.pdf') },
      content = function(file) {
        #ggsave(file = file, plot = plotInput(), device = "pdf")
        pdf(file)
        plotInput()
        dev.off()
      }
    )
    
    
  }
)
