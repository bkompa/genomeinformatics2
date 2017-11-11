ParseHMMScanFile <- function(filePath){
  person <- unlist(strsplit(getwd(), .Platform$file.sep))[2] 
  if(person!='scott' & person!='Mack') {
    .libPaths('/home/zmx21/R/x86_64-pc-linux-gnu-library/3.1')
    .libPaths('/home/asn32/R/x86_64-pc-linux-gnu-library/3.1')
  }
  
  data <- readLines(con=file(filePath),skipNul = T)
  
  queries <- data[grep(pattern = 'Query:',data)]
  header<- grep(pattern="--- full sequence ---",data)+3
  footer<- grep(pattern="Domain annotation for each model:",data)-3

  result <- data.frame(fullSeqEVal=c(),fullSeqScore=c(),fullSeqBias=c(),domEVal=c(),domScore=c(),domBias=c(),
                       domExp=c(),domN=c(),domModel=c(),domDescrip=c())
  
  for (i in 1:length(header)) {
    print(i)
    currentLines<-data[header[i]:footer[i]]
    currentLines<-currentLines[!grepl(pattern='  ------ inclusion threshold ------',currentLines)]
    splitLines <- mapply(strsplit,currentLines,' ')
    
    #If no hits returned
    if (length(splitLines[[1]]) == 0){
      next
    }
    
    fullSeqEVal <- rep(NA,length(splitLines))
    fullSeqScore <- rep(NA,length(splitLines))
    fullSeqBias <- rep(NA,length(splitLines))
    domEVal <- rep(NA,length(splitLines))
    domScore <- rep(NA,length(splitLines))
    domBias <- rep(NA,length(splitLines))
    domExp <- rep(NA,length(splitLines))
    domN <- rep(NA,length(splitLines))
    domModel <- vector(mode='character',length=length(splitLines))
    domDescrip <- vector(mode='character',length=length(splitLines))

    for (j in 1:length(splitLines)){
      currentLine <- splitLines[[j]]
      
      currentLine <- currentLine[!mapply(grepl,currentLine,"",fixed=T)]   #Pick out all elements that's not empty
      
      fullSeqEVal[j] <- as.numeric(currentLine[1]);fullSeqScore[j] <- as.numeric(currentLine[2]);
      fullSeqBias[j] <- as.numeric(currentLine[3]);domEVal[j] <- as.numeric(currentLine[4])
      domScore[j] <- as.numeric(currentLine[5]);domBias[j] <- as.numeric(currentLine[6]);
      domExp[j] <- as.numeric(currentLine[7]); domN[j] <- as.numeric(currentLine[8])
      domModel[j] <- currentLine[9]
      domDescrip[j] <- paste(currentLine[10:length(currentLine)],collapse = " ")
    }
    currentQuery <- unlist(strsplit(queries[i]," "))
    currentQuery <- currentQuery[!mapply(grepl,currentQuery,"",fixed=T)]
    qLength <- as.numeric(gsub(x=gsub("\\[L=","",currentQuery[3]),pattern = "\\]",replacement = ""))
    
    currentResult <- data.frame(qID = currentQuery[2],qLength = qLength,fullSeqEVal=fullSeqEVal,fullSeqScore=fullSeqScore,fullSeqBias=fullSeqBias,
                                domEVal=domEVal,domScore=domScore,domBias=domBias,
                                domExp=domExp,domN=domN,domModel=domModel,domDescrip=domDescrip)
    result <- rbind(result,currentResult)
    
  }
  eval(parse(text=paste0(filePath,'=result'))) #Reassign variable
  str = paste0("save(",filePath,",file=\'",filePath,".dat\')")   #Save file
  eval(parse(text=str))
}

ImportGODb <- function(filePath){
  db <- readLines(con=file(filePath),skipNul = T)
  db <- db[grepl('Pfam:',db)]  #Remove headers
  
  GoIds <- as.character(mapply(sub,pattern=' GO:',replacement='',x=as.character(sapply(db,function(x) unlist(strsplit(x,';'))[2]))))
  GoDescrip <- as.character(mapply(sub,pattern=' GO:',replacement='',x=as.character(sapply(db,function(x) unlist(strsplit(unlist(strsplit(x[1],'>'))[2],';'))[1]))))
  PfamId <- as.character(mapply(sub,pattern='Pfam:',replacement='',x=as.character(sapply(db,function(x) unlist(strsplit(x,' '))[1]))))
  return(data.frame(PfamId = PfamId,GoIDs=GoIds,GoDescrip=GoDescrip))
}

GetGoIDFromPfam <- function(PfamID,GODb){
  
}
