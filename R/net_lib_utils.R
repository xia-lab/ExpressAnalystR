##################################################
## R script for NetworkAnalyst
## Description: Functions to load various libraries for functional enrichment analysis during network visualization
## Authors:
## G. Zhou (guangyan.zhou@mail.mcgill.ca)
## J. Xia, jeff.xia@mcgill.ca
###################################################

# note, last two par only for STRING database
QueryPpiSQLite <- function(sqlite.path, table.nm, q.vec, requireExp, min.score){
  require('RSQLite')
  db.path <- paste(sqlite.path, "ppi.sqlite", sep="");
  ppi.db <- .connect.sqlite(db.path);
  query <- paste(shQuote(q.vec),collapse=",");
  
  if(grepl("string$", table.nm)){
    if(requireExp){
      statement <- paste("SELECT * FROM ",table.nm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))  AND combined_score >=", min.score, " AND experimental > 0", sep="");
    }else{
      statement <- paste("SELECT * FROM ",table.nm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))  AND combined_score >=", min.score, sep="");
    }
  }else{
    statement <- paste("SELECT * FROM ",table.nm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))", sep="");
  }
  ppi.res <- .query.sqlite(ppi.db, statement);
  
  # remove dupliated edges
  ppi.res <- ppi.res[!duplicated(ppi.res$row_id),]
  return(ppi.res);  
}

#for signaling pathway as well
QueryOtherPpiSQLite <- function(sqlite.path, table.nm, q.vec, data.org){
  require('RSQLite')
    if(grepl("signal$", table.nm)){
        db.path <- paste(sqlite.path, "signor.sqlite", sep=""); 
        table.nm <- data.org
    }else{
        db.path <- paste(sqlite.path, "cross_species_ppi.sqlite", sep=""); 
    }
  con <- .connect.sqlite(db.path);
  query <- paste(shQuote(q.vec),collapse=",");
  statement <- paste('SELECT * FROM "',table.nm, '" WHERE ((id1 IN (', query, ')) OR (id2 IN (', query, ')))', sep='');
  ppi.res <- .query.sqlite(con, statement);
  return(ppi.res);  
}

# table name is org code, id.type is column name
QueryMirSQLite <- function(sqlite.path, org, id.type, q.vec, db.nm){
  require('RSQLite');
  db.path <- paste(sqlite.path, "mir2gene.sqlite", sep=""); 
  con <- .connect.sqlite(db.path);
  query <- paste (shQuote(q.vec),collapse=",");
  #statement <- paste("SELECT * FROM ", org, " WHERE ",id.type," IN (",query,")", " AND ", db.nm," == 1 ", sep="");
  statement <- paste("SELECT * FROM ", org, " WHERE ",id.type," IN (",query,")", " AND ", db.nm," == 1 ", sep="");
  return(.query.sqlite(con, statement));
}

# table name is org code, id.type is column name
QueryDrugSQLite <- function(sqlite.path, q.vec){
  require('RSQLite');
  db.path <- paste(sqlite.path, "drug.sqlite", sep="");
  con <- .connect.sqlite(db.path);
  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM human WHERE upid IN (",query,")", sep="");
  return(.query.sqlite(con, statement));
}

QueryDiseaseSQLite <- function(sqlite.path, q.vec){
  require('RSQLite');
  db.path <- paste(sqlite.path, "disease.sqlite", sep="");
  con <- .connect.sqlite(db.path);
  query <- paste(shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM human WHERE entrez IN (",query,")", sep="");
  return(.query.sqlite(con, statement));
}

QueryTfmirSQLite <- function(sqlite.path, q.vec, data.org){
  require('RSQLite');
  db.path <- paste(sqlite.path, "tfmirgene.sqlite", sep="");
  con <- .connect.sqlite(db.path);
  query <- paste(shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", data.org, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))" , sep="");
  return(.query.sqlite(con, statement));
}

QueryDiffNetSQLite <- function(sqlite.path, q.vec){
  require('RSQLite');
  db.path <- paste(sqlite.path, "tissuePPI.sqlite", sep="");
  con <- .connect.sqlite(db.path);
  table.nm <- diffNetName;
  query <- paste (shQuote(q.vec),collapse=",");
  topPct <- 1-diffPct;
  botstatement <- paste("SELECT * FROM ",table.nm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))  AND rank <=", diffPct ,sep="");
  topstatement <- paste("SELECT * FROM ",table.nm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))  AND rank >=", topPct ,sep="");
  
  dic1 <- .query.sqlite(con, botstatement, FALSE);# no close db connection
  dic2 <- .query.sqlite(con, topstatement);
  dic <- rbind(dic1, dic2);
  return(dic);
}

QueryCellCoexSQLite <- function(sqlite.path, q.vec, data.org){
  require('RSQLite');
  db.path <- paste(sqlite.path, data.org,"_immune.sqlite", sep="");
  con <- .connect.sqlite(db.path);
  tblNm <- paste0(data.org,"_",cellCoexNumber);
  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", tblNm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))" , sep="");
  return(.query.sqlite(con,statement));
}

QueryTissueCoexSQLite <- function(sqlite.path, q.vec, data.org){
  require('RSQLite');
  db.path <- paste(sqlite.path, "tissueCoex.sqlite", sep="");
  con <- .connect.sqlite(db.path);
  tblNm <- paste0(data.org,"_",cellCoexNumber);
  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", tblNm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))" , sep="");
  return(.query.sqlite(con, statement));
}

QueryChemSQLite<- function(sqlite.path, org, q.vec){
  require('RSQLite');
  db.path <- paste(sqlite.path, "chem.sqlite", sep="");
  con <- .connect.sqlite(db.path);
  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", org, " WHERE entrez IN (",query,")", sep="");
  return(.query.sqlite(con, statement));
}

QueryTFSQLite<- function(sqlite.path, table.nm, q.vec){
  require('RSQLite');
  db.path <- paste(sqlite.path, "tf2gene.sqlite", sep="");
  con <- .connect.sqlite(db.path);
  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", table.nm, " WHERE entrez IN (",query,")", sep="");
  return(.query.sqlite(con, statement));
}

doPpiIDMapping <- function(sqlite.path, q.vec, data.org="entrez_swissprot"){
  paramSet <- readSet(paramSet, "paramSet");
  data.org <- paramSet$data.org;
  invertAcc <-F;
  idType <- ""
  if(data.org == "sce"){
    if(net.type == "string"){ # only for yeast
      idType <- "entrez";
    }else{
      idType <- "entrez_uniprot";
    }
  }else if(net.type %in% c("irefinx", "rolland")){
    idType <- "entrez_uniprot";
  }else{
    if(data.org == "pae" & net.type == "interactome"){
      idType <- "entrez_string";
    }else if(data.org %in% c("tbr", "cel", "pfa")){
      idType <- "entrez_string";
    }else{
      idType <- "entrez";
      invertAcc<- T;
    }
  }
  db.map <-  queryGeneDB(idType, data.org);
  hit.inx <- match(q.vec, db.map[, "gene_id"]);
  ppi.mat <- db.map[hit.inx, ];

if(idType == "entrez"){
  return(ppi.mat$gene_id);
  }else{
  return(ppi.mat$accession);
}
}

doEntrez2UniprotMapping<-function(entrez.vec, paramSet){
  data.org <- paramSet[["data.org"]];
 
  db.map <-  queryGeneDB("entrez_uniprot", data.org);
  hit.inx <- match(entrez.vec, db.map[, "gene_id"]);
  entrezs <- db.map[hit.inx, "accession"];
  mode(entrezs) <- "character";
  na.inx <- is.na(entrezs);
  entrezs[na.inx] <- entrez.vec[na.inx];
  return(entrezs);
}

QueryMetSQLiteNet <- function(sqlite.path, table.nm, q.vec, inv){
    require('RSQLite');
    db.path <- paste0(sqlite.path, "met.sqlite")
    con <- .connect.sqlite(db.path);
    query <- paste (shQuote(q.vec),collapse=",");
    
    if(inv == "inverse"){
        statement <- paste("SELECT * FROM ", table.nm, " WHERE entrez IN (",query,")", sep="");
    }else{
        statement <- paste("SELECT * FROM ", table.nm, " WHERE kegg IN (",query,")", sep="");
    } 

    lnctable <- dbSendQuery(con, statement);
    lnc.dic <- fetch(lnctable, n=-1); # get all records
    dbDisconnect(con);
    return(lnc.dic);
}

##########################################
############# private utility methods #### 
##########################################

.query.sqlite <- function(db.con, statement, offline=TRUE){
  rs <- dbSendQuery(db.con, statement);
  res <- fetch(rs, n=-1); # get all records
  dbClearResult(rs);
  if(offline){
    dbDisconnect(db.con);
  }
  cleanMem();
  return(res);
}

.connect.sqlite <- function(db.path){
  if(!PrepareSqliteDB(db.path, paramSet$on.public.web)){
    stop("Sqlite database is missing, please check your internet connection!");
  }
  return(dbConnect(SQLite(), db.path)); 
}

PrepareSqliteDB <- function(sqlite_Path, onweb = TRUE) {
  if(onweb) {return(TRUE)};
  if(file.exists(sqlite_Path)) {return(TRUE)};

  dbNM <- basename(sqlite_Path);
  DonwloadLink <- paste0("https://www.xialab.ca/resources/sqlite/", dbNM);
  download.file(DonwloadLink, sqlite_Path);
  return(TRUE)
}

queryGeneDB <- function(table.nm, data.org){
  paramSet <- readSet(paramSet, "paramSet");  
  if(length(table.nm) == 0){
    table.nm <- "";
  }

  if(table.nm == "custom" || data.org == "custom"){
    db.map <- qs::qread("anot_table.qs");
  }else{
    require('RSQLite');
    db.path <- paste(paramSet$sqlite.path, data.org, "_genes.sqlite", sep="")

    if(!PrepareSqliteDB(db.path, paramSet$on.public.web)){
      stop("Sqlite database is missing, please check your internet connection!");
    }
    conv.db <- dbConnect(SQLite(), db.path); 
    tbls <- dbListTables(conv.db)
    if(!table.nm %in% tbls){
        return(0);
    }
    db.map <- dbReadTable(conv.db, table.nm);
    dbDisconnect(conv.db); cleanMem();
  }
  return(db.map)
}