
WriteSampleTable <- function(ftpDir, userDir, sampleTable){
  str = gsub("\\[|\\]", "", sampleTable);
  str = gsub("FastqModel\\{|\\}", "", str);
  
  tmp = unlist(strsplit(str, ", "));
  nSamples = length(tmp) / 9;
  paired = sum(tmp %in% "paired=true") > 0;

  # OPTIMIZED: Vectorized operations to eliminate O(n²) memory allocations
  # Pre-compute indices once
  indices <- (seq_len(nSamples) - 1) * 9

  # Vectorized string operations (no loop!)
  group <- gsub("group=", "", tmp[indices + 1])
  baseNamesSim <- gsub("pseudoName=", "", tmp[indices + 2])
  baseNames <- file.path(userDir, baseNamesSim)
  fReadFiles <- file.path(ftpDir, gsub("forwardReadFile=", "", tmp[indices + 3]))

  if(paired){
    rReadFiles <- file.path(ftpDir, gsub("reverseReadFile=", "", tmp[indices + 4]))
  } else {
    rReadFiles <- c()
  }
  
  if(paired){
    samTable = data.frame("sampleName" = baseNames,
                          "forwardReads" = fReadFiles,
                          "reverseReads" = rReadFiles);
  } else {
    samTable = data.frame("sampleName" = baseNames,
                          "forwardReads" = fReadFiles);
  }
  
  samTable$group <- group;
  sampGroup = data.frame("sampleName" = baseNamesSim,
                         "sampleGroup" = group);
  
  write.table(samTable,
              file.path(userDir, "sampleTable.txt"),
              quote = F,
              row.names = F,
              sep = "\t",
              col.names = F);
  
  write.table(sampGroup,
              file.path(userDir, "sampleGroup.txt"),
              quote = F,
              row.names = F,
              sep = "\t",
              col.names = F);
  
  return (1);
  
}

WriteJobId <- function(userDir, email, jobId){
  write.table(jobId,
              file.path(userDir, email, "jobId.txt"),
              quote = F,
              row.names = F,
              sep = "\t",
              col.names = F);
  
  #jobId <- read.table(file.path(userDir, email, "jobId.txt"));
  
}

SubmitJobS2f <- function(userDir, email, database, mismatch, minlength, minscore, des, outputmappedreads, shellscriptDir){
  
  isSlurm = FALSE;
  isPengLocal = F;
  if(file.exists("/home/zgy/NetBeansProjects/")){
    seq2funPath = "/data/glassfish/seq_software/seq2fun_source/Seq2Fun/bin/seq2fun";
    databasePath = "/data/glassfish/seq_software/seq2fun_database";
    isSlurm = FALSE;
  } else if(grepl("glassfish", userDir)){
    seq2funPath = "/data/glassfish/seq_software/seq2fun_source/Seq2Fun/bin/seq2fun";
    databasePath = "/data/glassfish/seq_software/seq2fun_database";
    isSlurm = TRUE;
  } else if(grepl("peng", userDir)) {
    seq2funPath = "/home/peng/expressanalyst_test/software/Seq2Fun/bin/seq2fun";
    databasePath = "/home/peng/expressanalyst_test/seq2fun_database";
    isSlurm = FALSE;
  } else if(grepl("jeffxia/Dropbox", userDir)) { # jeff local
    seq2funPath = "/Users/jeffxia/Dropbox/resources/eoa/software/Seq2Fun/bin/seq2fun";
    databasePath = "/Users/jeffxia/Dropbox/resources/eoa/seq2fun_database";
    isSlurm = FALSE;
  } else if(grepl("ewald", userDir)) { # jess local
    seq2funPath = "/Users/jessicaewald/NetbeansProjects/Seq2Fun/bin/seq2fun";
    databasePath = "/Users/jessicaewald/eoa/seq2fun_database/";
    isSlurm = F;
  } else {
    
  }
  
  if(isSlurm){
    # OPTIMIZED: Single paste() call instead of repeated paste0() to eliminate O(n²) string copying
    conf_inf <- paste(
      "#!/bin/bash",
      "#",
      "#SBATCH --job-name=Seq2fun_Processing",
      "#",
      "#SBATCH --ntasks=1",
      "#SBATCH --time=14400:00",
      "#SBATCH --mem-per-cpu=5000",
      "#SBATCH --cpus-per-task=4",
      paste0("#SBATCH --mail-user=", email),
      "#SBATCH --mail-type=BEGIN",
      "#SBATCH --mail-type=END",
      "#SBATCH --mail-type=FAIL",
      "#SBATCH --mail-type=REQUEUE",
      "#SBATCH --mail-type=ALL",
      paste0("#SBATCH --output=", file.path(userDir, "slurm_logout.txt")),
      sep = "\n"
    )

    if(file.exists("/docker_marker")){
      conf_inf <- paste(
      "#!/bin/bash",
      "#",
      "#SBATCH --job-name=Seq2fun_Processing",
      "#",
      "#SBATCH --ntasks=1",
      "#SBATCH --time=28800:00",
      paste0("#SBATCH --output=", file.path(userDir, "slurm_logout.txt")),
      sep = "\n"
    )
    }
    
    str <- paste(seq2funPath,
                 "--longlog --sampletable",
                 file.path(userDir, "sampleTable.txt"),
                 "--tfmi",
                 file.path(databasePath, paste0(database, "_v*.fmi")),
                 "--genemap",
                 file.path(databasePath, paste0(database, "_annotation_v*.txt")),
                 "--mismatch", mismatch,
                 "--minlength", minlength,
                 "--minscore", minscore,
                 "-w 4 --profiling -V;",
                 sep = " ");
    
    strcmd <- paste0("seq2fun: ",
                     "\ndatabase: ", database,
                     "\nmaximum number of mismatches: ", mismatch,
                     "\nminimum matching length: ", minlength,
                     "\nminimum matching blosum62 score: ", minscore,
                     "\nproject_description: ", des);
   
    
    str_R <- paste("Rscript --vanilla",
                   file.path(shellscriptDir, "run_process.R"),
                   userDir,
                   shellscriptDir,
                   sep = " ");
    
    str_checMysql1 <- paste("if ! [[ -x ",
                            file.path(userDir, "renameDirFromMysql.sh"),
                            " ]]; then chmod +x ",
                            file.path(userDir, "renameDirFromMysql.sh"),
                            "; fi;",
                            sep = " ");
    
    str_runMysql1 <- paste("bash ",
                           file.path(userDir, "renameDirFromMysql.sh"),
                           ";",
                           sep = " ");
    
    if(isPengLocal){
      str_checMysql2 <- paste("if ! [[ -x ",
                              file.path(userDir, "updateMysql.sh"),
                              " ]]; then chmod +x ",
                              file.path(userDir, "updateMysql.sh"),
                              "; fi;",
                              sep = " ");
      
      str_runMysql2 <- paste("bash ",
                             file.path(userDir, "updateMysql.sh"),
                             "; rm -rf ",
                             file.path(userDir, "updateMysql.sh"),
                             " 2>&1 | tee -a ",
                             file.path(userDir, "slurm_logout.txt"),
                             sep = " ");
    }
    
    sink(file.path(userDir, "ExecuteRawSeq.sh"));
    cat(conf_inf, "\n\n");
    cat(str, "\n\n");
    cat(str_R, "\n\n");
    if(isPengLocal){
      cat(str_checMysql2, "\n\n");
      cat(str_runMysql2, "\n\n");
    }
    
    # cat(str_runRemove, "\n\n");
    # cat(str_zipDir, "\n\n");
    cat(str_checMysql1, "\n\n");
    cat(str_runMysql1, "\n\n");
    sink();
    
    sink(file.path(userDir, "analysis_parameters.txt"));
    cat(strcmd, "\n");
    sink();
    
  } else {

    str <- paste(seq2funPath,
                 "--longlog --sampletable",
                 file.path(userDir, "sampleTable.txt"),
                 "--tfmi",
                 file.path(databasePath, paste0(database, "_v*.fmi")),
                 "--genemap",
                 file.path(databasePath, paste0(database, "_annotation_v*.txt")),
                 "--mismatch", mismatch,
                 "--minlength", minlength,
                 "--minscore", minscore,
                 "-w 4 --profiling -V;",
                 sep = " ");
    
    strcmd <- paste0("seq2fun: ",
                     "\ndatabase: ", database,
                     "\nmaximum number of mismatches: ", mismatch,
                     "\nminimum matching length: ", minlength,
                     "\nminimum matching blosum62 score: ", minscore,
                     "\nproject_description: ", des);
        
    str_R <- paste("Rscript --vanilla",
                   file.path(shellscriptDir, "run_process.R"),
                   userDir,
                   shellscriptDir,
                   sep = " ");
    
    str_checMysql1 <- paste("if ! [[ -x ",
                            file.path(userDir, "renameDirFromMysql.sh"),
                            " ]]; then chmod +x ",
                            file.path(userDir, "renameDirFromMysql.sh"),
                            "; fi;",
                            sep = " ");
    
    str_runMysql1 <- paste("bash ",
                           file.path(userDir, "renameDirFromMysql.sh"),
                           ";",
                           sep = " ");
    
    sink(file.path(userDir, "ExecuteRawSeq.sh"));
    cat(str, "\n\n");
    cat(str_R, "\n\n");
    
    cat(str_checMysql1, "\n\n");
    cat(str_runMysql1, "\n\n");
    sink();
    
    sink(file.path(userDir, "analysis_parameters.txt"));
    cat(strcmd, "\n");
    sink();
    
    str_sh <- "#!/bin/sh\n";
    
    str_sh <- paste0(str_sh,
                     "chmod +x ",
                     file.path(userDir, "ExecuteRawSeq.sh; "),
                     "nohup bash ",
                     file.path(userDir, "ExecuteRawSeq.sh > "),
                     file.path(userDir, "slurm_logout.txt 2>&1 &\n"),
                     "echo $! > ",
                     file.path(userDir, "jobid_pid.txt;"));
    
    sink(file.path(userDir, "launcher.sh"));
    cat(str_sh, "\n");
    sink();
    
  }
  return(1);
}


submitMysqlJobRename <- function(userDir, jobID, aligner){
  pDir = dirname(userDir);
  
  str_mysql1 = "";
  if(aligner == "Seq2Fun"){
    str_mysql1 <- paste(
      "echo Starting zip files $(date);",
      "zip -r -j",
      file.path(userDir, "Download.zip"),
      file.path(userDir, "S2fid_*.txt"),
      file.path(userDir, "Seq2Fun_summary_all_samples.html"),
      file.path(userDir, "*.png"),
      file.path(userDir, "metaTable.txt"),
      file.path(userDir, "analysis_parameters.txt"),
      "-x",
      file.path(userDir, "*.fastq.gz;"),
      "echo Finish zip files $(date);",
      sep = " ");
  } else if(aligner == "DePlexer") {
    deplexJobOutput <- file.path(pDir, paste0("job_", jobID), "output")
    str_mysql1 <- paste(
      "echo Starting zip files $(date);",
      "zip -r",
      file.path(userDir, "Download.zip"),
      file.path(userDir, "output"),
      deplexJobOutput,
      file.path(userDir, "*_deplexer_out"),
      file.path(userDir, "*.html"),
      file.path(userDir, "analysis_parameters.txt"),
      "-x",
      file.path(userDir, "*.fastq.gz;"),
      "echo Finish zip files $(date);",
      sep = " ")
  } else {
    str_mysql1 <- paste(
      "echo Starting zip files $(date);",
      "zip -r -j",
      file.path(userDir, "Download.zip"),
      file.path(userDir, "All_samples_*_txi_counts.txt"),
      file.path(userDir, "*.png"),
      file.path(userDir, "metaTable.txt"),
      file.path(userDir, "analysis_parameters.txt"),
      "-x",
      file.path(userDir, "*.fastq.gz"),
      file.path(userDir, "All_samples_*txi_abundance.txt"),
      file.path(userDir, "All_samples_*txi_length.txt;"),
      "echo Finish zip files $(date);",
      sep = " ")
  }

  str_mysql2 <- paste("cp -r",
                     file.path(userDir, "*"),
                     file.path(pDir, paste0("job_", jobID)),
                     sep = " ")
  
  sink(file.path(userDir, "renameDirFromMysql.sh"));
  cat(str_mysql1, "\n");
  cat(str_mysql2, "\n");
  sink();
}

submitMysqlJobRenamePro <- function(userDir, jobID, aligner){
  pDir = dirname(userDir);
  
  str_mysql1 = "";
  if(aligner == "Seq2Fun"){
    str_mysql1 <- paste(
      "echo Starting zip files $(date);",
      "zip -r -j",
      file.path(userDir, "Download.zip"),
      file.path(userDir, "S2fid_*.txt"),
      file.path(userDir, "Seq2Fun_summary_all_samples.html"),
      file.path(userDir, "*.png"),
      file.path(userDir, "metaTable.txt"),
      file.path(userDir, "analysis_parameters.txt"),
      "-x",
      file.path(userDir, "*.fastq.gz;"),
      "echo Finish zip files $(date);",
      sep = " ");
  } else if(aligner == "DePlexer") {
    deplexJobOutput <- file.path(pDir, paste0("job_", jobID), "output")
    str_mysql1 <- paste(
      "echo Starting zip files $(date);",
      "zip -r",
      file.path(userDir, "Download.zip"),
      file.path(userDir, "output"),
      deplexJobOutput,
      file.path(userDir, "*_deplexer_out"),
      file.path(userDir, "*.html"),
      file.path(userDir, "analysis_parameters.txt"),
      "-x",
      file.path(userDir, "*.fastq.gz;"),
      "echo Finish zip files $(date);",
      sep = " ")
  } else {
    str_mysql1 <- paste(
      "echo Starting zip files $(date);",
      "zip -r -j",
      file.path(userDir, "Download.zip"),
      file.path(userDir, "All_samples_*_txi_counts.txt"),
      file.path(userDir, "*.png"),
      file.path(userDir, "metaTable.txt"),
      file.path(userDir, "analysis_parameters.txt"),
      "-x",
      file.path(userDir, "*.fastq.gz"),
      file.path(userDir, "All_samples_*txi_abundance.txt"),
      file.path(userDir, "All_samples_*txi_length.txt;"),
      "echo Finish zip files $(date);",
      sep = " ")
  }

  sink(file.path(userDir, "renameDirFromMysql.sh"));
  cat(str_mysql1, "\n");
  sink();
}

submitMysqlJobStatus <- function(userDir, jobID, numSamples){
  str_mysql <- paste("mysql -u proftpd -pseq2fun ftp -e ",
                     "\"UPDATE jobs SET jobStatus = 'COMPLETED', numFinSamples =",
                     numSamples,
                     "WHERE jobID = ",
                     jobID,
                     "\";",
                     sep = " ")
  
  sink(file.path(userDir, "updateMysql.sh"));
  cat(str_mysql, "\n");
  sink();
}

SubmitJobSlm <- function(userDir, email, database, des, readEnds, shellscriptDir, minScore){
  isPengLocal = FALSE;
   if(file.exists("/home/zgy/NetBeansProjects/")){
    fastpPath = "/data/glassfish/projects/expressanalyst/fastp_source/fastp";
    salmonPath = "/data/glassfish/projects/expressanalyst/salmon_source/salmon-latest_linux_x86_64/bin/salmon";
    databasePath = "/data/glassfish/projects/expressanalyst/salmon_database";
    isPengLocal = TRUE;
  } else if(grepl("glassfish", userDir)){
    fastpPath = "/data/glassfish/projects/expressanalyst/fastp_source/fastp";
    salmonPath = "/data/glassfish/projects/expressanalyst/salmon_source/salmon-latest_linux_x86_64/bin/salmon";
    databasePath = "/data/glassfish/projects/expressanalyst/salmon_database";
    isPengLocal = FALSE;
  } else {
    fastpPath = "/home/peng/software/fastp/fastp";
    salmonPath = "/home/peng/software/salmon-latest_linux_x86_64/bin/salmon";
    databasePath = "/home/peng/expressanalyst_test/database/salmon";
    isPengLocal = TRUE;
  }
  
  fastpOutputR1 = file.path(userDir, "tmp_clean_R1.fastq.gz");
  if(readEnds == "pe"){
    fastpOutputR2 = file.path(userDir, "tmp_clean_R2.fastq.gz");
  }
  
  ## Prepare Configuration script for slurm running
  # OPTIMIZED: Single paste() call instead of repeated paste0() to eliminate O(n²) string copying
  conf_inf <- paste(
    "#!/bin/bash",
    "#",
    "#SBATCH --job-name=Salmon_Processing",
    "#",
    "#SBATCH --ntasks=1",
    "#SBATCH --time=14400:00",
    "#SBATCH --mem-per-cpu=5000",
    "#SBATCH --cpus-per-task=4",
    paste0("#SBATCH --mail-user=", email),
    "#SBATCH --mail-type=BEGIN",
    "#SBATCH --mail-type=END",
    "#SBATCH --mail-type=FAIL",
    "#SBATCH --mail-type=REQUEUE",
    "#SBATCH --mail-type=ALL",
    paste0("#SBATCH --output=", file.path(userDir, "slurm_logout.txt")),
    sep = "\n"
  )

    
    if(file.exists("/docker_marker")){
        conf_inf <- paste(
          "#!/bin/bash",
          "#",
          "#SBATCH --job-name=Salmon_Processing",
          "#",
          "#SBATCH --ntasks=1",
          "#SBATCH --time=28800:00",
          paste0("#SBATCH --output=", file.path(userDir, "slurm_logout.txt")),
          sep = "\n"
        )
    }

  
  if(readEnds == "pe"){
    str <- paste("cat",
                 file.path(userDir, "sampleTable.txt"),
                 "| while read fn ff rf grp; do echo processing sample ${fn};",
                 fastpPath, 
                 "--qualified_quality_phred", minScore,
                 "-i ${ff} -I ${rf} -o", fastpOutputR1, "-O", fastpOutputR2,
                 "-w 4 -h ${fn}_fastp.html -j ${fn}_fastp.json -c -y -V &&",
                 salmonPath, "quant -i",
                 file.path(databasePath, database, "salmon_index"),
                 "-l A -1", fastpOutputR1, "-2", fastpOutputR2, "-o ${fn}_salmon",
                 "-p 4 && ", "mv ${fn}_salmon/quant.sf ${fn}_salmon_quant.txt && ",
                 "mv ${fn}_salmon/logs/salmon_quant.log ${fn}_salmon_log.txt && ",
                 "rm -rf ${fn}_salmon && ",
                 "echo Fastp and Salmon finish sample ${fn}; done && rm",
                 file.path(userDir, "tmp_clean_*.gz"),
                 "&& echo Fastp and Salmon have processed all samples;",
                 sep = " ");
    
  } else {
    str <- paste("cat",
                 file.path(userDir, "sampleTable.txt"),
                 "| while read fn ff grp; do echo processing sample ${fn};",
                 fastpPath, 
                 "--qualified_quality_phred", minScore,
                 "-i ${ff} -o", fastpOutputR1,
                 "-w 4 -h ${fn}_fastp.html -j ${fn}_fastp.json -y &&",
                 salmonPath, "quant -i",
                 file.path(databasePath, database, "salmon_index"),
                 "-l A -r", fastpOutputR1, "-o ${fn}_salmon",
                 "-p 4 && ", "mv ${fn}_salmon/quant.sf ${fn}_salmon_quant.txt && ",
                 "mv ${fn}_salmon/logs/salmon_quant.log ${fn}_salmon_log.txt && ",
                 "rm -rf ${fn}_salmon && ",
                 "echo Fastp and Salmon finish sample ${fn}; done && rm",
                 file.path(userDir, "tmp_clean_*.gz"),
                 "&& echo Fastp and Salmon have processed all samples;",
                 sep = " ");
  }
  
strcmd <- paste0("fastp: ", 
                "\nmininum quality score: ", minScore,
                "\nsalmon: ", 
                "\ndatabase: ", list.files(file.path(databasePath, database), pattern = ".gz$"),
                "\nproject_description: ", des);

  str_R <- paste("echo Starting process tables and figures in R $(date);",
                 "Rscript --vanilla",
                 file.path(shellscriptDir, "run_process_salmon.R"),
                 userDir,
                 "; echo Finish process tables and figures in R $(date);",
                 sep = " ");

  str_checMysql1 <- paste("if ! [[ -x ",
                          file.path(userDir, "renameDirFromMysql.sh"),
                          " ]]; then chmod +x ",
                          file.path(userDir, "renameDirFromMysql.sh"),
                          "; fi;",
                          sep = " ");
  
  str_runMysql1 <- paste("bash ",
                         file.path(userDir, "renameDirFromMysql.sh"),
                         ";",
                         sep = " ");
  
  # str_runRemove <- paste("mkdir -p",
  #                        file.path(dirname(userDir), "results;"),
  #                        "cp -r",
  #                        file.path(dirname(userDir), "All_sample*counts.txt"),
  #                        file.path(dirname(userDir), "*png"),
  #                        file.path(dirname(userDir), "metaTable.txt"),
  #                        file.path(dirname(userDir), "results;"),
  #                        sep = " ");
  # 
  # str_zipDir <- paste("echo Starting zip files $(date);",
  #                     "zip -r -j",
  #                     file.path(userDir, "Download.zip"),
  #                     file.path(userDir, "results"),
  #                     "; echo Finish zip files $(date);",
  #                     sep = " ");
  
  if(isPengLocal){
    str_checMysql2 <- paste("if ! [[ -x ",
                            file.path(userDir, "updateMysql.sh"),
                            " ]]; then chmod +x ",
                            file.path(userDir, "updateMysql.sh"),
                            "; fi;",
                            sep = " ");
    
    str_runMysql2 <- paste("bash ",
                           file.path(userDir, "updateMysql.sh"),
                           "; rm -rf ",
                           file.path(userDir, "updateMysql.sh"),
                           file.path(userDir, "tmp_clean_R1.fastq.gz"),
                           file.path(userDir, "*_salmon"),
                           ";",
                           sep = " ");
  }
  
  sink(file.path(userDir, "ExecuteRawSeq.sh"));
  cat(conf_inf, "\n\n");
  cat(str, "\n\n");
  cat(str_R, "\n\n");
  if(isPengLocal){
    cat(str_checMysql2, "\n\n");
    cat(str_runMysql2, "\n\n");
  }
  # cat(str_runRemove, "\n\n");
  # cat(str_zipDir, "\n\n");
  cat(str_checMysql1, "\n\n");
  cat(str_runMysql1, "\n\n");
  sink();
    sink(file.path(userDir, "analysis_parameters.txt"));
    cat(strcmd, "\n");
    sink();
  return(1);
}


SubmitJobKls <- function(userDir, email, database, des, readEnds, shellscriptDir, minScore, aveFragLen, stdFragLen){
  userDir <<- userDir;
  email <<- email;
  database <<- database;
  des <<- des;
  readEnds <<- readEnds;
  shellscriptDir <<- shellscriptDir;
  minScore <<- minScore;
  aveFragLen <<- aveFragLen;
  stdFragLen <<- stdFragLen;
  
  isSlurm = FALSE;
  isPengLocal = F;
  if(file.exists("/home/zgy/NetBeansProjects/")){
    fastpPath = "/data/glassfish/seq_software/fastp_source/fastp";
    kallistoPath = "/data/glassfish/seq_software/kallisto_source/kallisto";
    databasePath = "/data/glassfish/seq_software/kallisto_database";
    isPengLocal = FALSE;
    isSlurm = F;
  }else if(grepl("glassfish", userDir)){
    fastpPath = "/data/glassfish/seq_software/fastp_source/fastp";
    kallistoPath = "/data/glassfish/seq_software/kallisto_source/kallisto";
    databasePath = "/data/glassfish/seq_software/kallisto_database";
    isPengLocal = FALSE;
    isSlurm = T;
  }else if(grepl("jeffxia/Dropbox", userDir)) { # jeff local
    fastpPath = "/Users/xia/Dropbox/resources/eoa/software/fastp_source/fastp";
    kallistoPath = "/Users/xia/Dropbox/resources/eoa/software/kallisto_source/kallisto";
    databasePath = "/Users/xia/Dropbox/resources/eoa/kallisto_database/";
    isSlurm = F;
  } else {
    fastpPath = "/home/peng/software/fastp/fastp";
    kallistoPath = "/home/peng/software/kallisto/kallisto";
    databasePath = "/home/peng/expressanalyst_test/database/kallisto";
    isPengLocal = TRUE;
  }
  
  fastpOutputR1 = file.path(userDir, "tmp_clean_R1.fastq.gz");
  if(readEnds == "pe"){
    fastpOutputR2 = file.path(userDir, "tmp_clean_R2.fastq.gz");
  }
  if(isSlurm){
    ## Prepare Configuration script for slurm running
    # OPTIMIZED: Single paste() call instead of repeated paste0() to eliminate O(n²) string copying
    conf_inf <- paste(
      "#!/bin/bash",
      "#",
      "#SBATCH --job-name=Salmon_Processing",
      "#",
      "#SBATCH --ntasks=1",
      "#SBATCH --time=14400:00",
      "#SBATCH --mem-per-cpu=5000",
      "#SBATCH --cpus-per-task=4",
      paste0("#SBATCH --mail-user=", email),
      "#SBATCH --mail-type=BEGIN",
      "#SBATCH --mail-type=END",
      "#SBATCH --mail-type=FAIL",
      "#SBATCH --mail-type=REQUEUE",
      "#SBATCH --mail-type=ALL",
      paste0("#SBATCH --output=", file.path(userDir, "slurm_logout.txt")),
      sep = "\n"
    )
    
    if(file.exists("/docker_marker")){
    conf_inf <- paste(
      "#!/bin/bash",
      "#",
      "#SBATCH --job-name=Salmon_Processing",
      "#",
      "#SBATCH --ntasks=1",
      "#SBATCH --time=28800:00",
      paste0("#SBATCH --output=", file.path(userDir, "slurm_logout.txt")),
      sep = "\n"
    )
    }

    if(readEnds == "pe"){
      str <- paste("cat",
                   file.path(userDir, "sampleTable.txt"),
                   "| while read fn ff rf grp; do echo processing sample ${fn};",
                   fastpPath, 
                   "--qualified_quality_phred", minScore,
                   "-i ${ff} -I ${rf} -o", fastpOutputR1, "-O", fastpOutputR2,
                   "-w 4 -h ${fn}_fastp.html -j ${fn}_fastp.json -c -y -V &&",
                   kallistoPath, "quant -i",
                   file.path(databasePath, database, "index/transcripts.idx"),
                   "-o ${fn}_kallisto -b 100 -t 4",
                   fastpOutputR1, fastpOutputR2, 
                   " && mv ${fn}_kallisto/abundance.tsv ${fn}_kallisto_quant.tsv && ", 
                   " mv ${fn}_kallisto/run_info.json ${fn}_kallisto_log.json && ",
                   "rm -rf ${fn}_kallisto && ",
                   "echo Fastp and Kallisto finish sample ${fn}; done && rm",
                   file.path(userDir, "tmp_clean_*.gz"),
                   "&& echo Fastp and Kallisto have processed all samples;",
                   sep = " ");
      
    } else {
      str <- paste("cat",
                   file.path(userDir, "sampleTable.txt"),
                   "| while read fn ff grp; do echo processing sample ${fn};",
                   fastpPath, 
                   "--qualified_quality_phred", minScore,
                   "-i ${ff} -o", fastpOutputR1,
                   "-w 4 -h ${fn}_fastp.html -j ${fn}_fastp.json -y &&",
                   kallistoPath, "quant -i",
                   file.path(databasePath, database, "index/transcripts.idx"),
                   "-o ${fn}_kallisto -b 100 -t 4",
                   "--single -l", aveFragLen, 
                   "-s", stdFragLen,
                   fastpOutputR1,
                   " && mv ${fn}_kallisto/abundance.tsv ${fn}_kallisto_quant.tsv && ", 
                   " mv ${fn}_kallisto/run_info.json ${fn}_kallisto_log.json && ",
                   "rm -rf ${fn}_kallisto && ",
                   "echo Fastp and Kallisto finish sample ${fn}; done && rm",
                   file.path(userDir, "tmp_clean_*.gz"),
                   "&& echo Fastp and Kallisto have processed all samples;",
                   sep = " ");
    }
    
    strcmd <- paste0("fastp: ", 
                     "\nmininum quality score: ", minScore,
                     "\nkallisto: ", 
                     "\ndatabase: ", list.files(file.path(databasePath, database), pattern = ".gz$"),
                     "\nproject_description: ", des);
    
    str_R <- paste("echo Starting process tables and figures in R $(date);",
                   "Rscript --vanilla",
                   file.path(shellscriptDir, "run_process_kallisto.R"),
                   userDir,
                   "; echo Finish process tables and figures in R $(date);",
                   sep = " ");
    
    str_checMysql1 <- paste("if ! [[ -x ",
                            file.path(userDir, "renameDirFromMysql.sh"),
                            " ]]; then chmod +x ",
                            file.path(userDir, "renameDirFromMysql.sh"),
                            "; fi;",
                            sep = " ");
    
    str_runMysql1 <- paste("bash ",
                           file.path(userDir, "renameDirFromMysql.sh"),
                           ";",
                           sep = " ");
    
    if(isPengLocal){
      str_checMysql2 <- paste("if ! [[ -x ",
                              file.path(userDir, "updateMysql.sh"),
                              " ]]; then chmod +x ",
                              file.path(userDir, "updateMysql.sh"),
                              "; fi;",
                              sep = " ");
      
      str_runMysql2 <- paste("bash ",
                             file.path(userDir, "updateMysql.sh"),
                             "; rm -rf ",
                             file.path(userDir, "updateMysql.sh"),
                             file.path(userDir, "tmp_clean_R1.fastq.gz"),
                             file.path(userDir, "*_kallisto"),
                             ";",
                             sep = " ");
    }
    
    sink(file.path(userDir, "ExecuteRawSeq.sh"));
    cat(conf_inf, "\n\n");
    cat(str, "\n\n");
    cat(str_R, "\n\n");
    if(isPengLocal){
      cat(str_checMysql2, "\n\n");
      cat(str_runMysql2, "\n\n");
    }
    # cat(str_runRemove, "\n\n");
    # cat(str_zipDir, "\n\n");
    cat(str_checMysql1, "\n\n");
    cat(str_runMysql1, "\n\n");
    sink();
    
    sink(file.path(userDir, "analysis_parameters.txt"));
    cat(strcmd, "\n");
    sink();
  }else{
    if(readEnds == "pe"){
      str <- paste("cat",
                   file.path(userDir, "sampleTable.txt"),
                   "| while read fn ff rf grp; do echo processing sample ${fn};",
                   fastpPath, 
                   "--qualified_quality_phred", minScore,
                   "-i ${ff} -I ${rf} -o", fastpOutputR1, "-O", fastpOutputR2,
                   "-w 4 -h ${fn}_fastp.html -j ${fn}_fastp.json -c -y -V &&",
                   kallistoPath, "quant -i",
                   file.path(databasePath, database, "index/transcripts.idx"),
                   "-o ${fn}_kallisto -b 100 -t 4",
                   fastpOutputR1, fastpOutputR2, 
                   " && mv ${fn}_kallisto/abundance.tsv ${fn}_kallisto_quant.tsv && ", 
                   " mv ${fn}_kallisto/run_info.json ${fn}_kallisto_log.json && ",
                   "rm -rf ${fn}_kallisto && ",
                   "echo Fastp and Kallisto finish sample ${fn}; done && rm",
                   file.path(userDir, "tmp_clean_*.gz"),
                   "&& echo Fastp and Kallisto have processed all samples;",
                   sep = " ");
      
    } else {
      str <- paste("cat",
                   file.path(userDir, "sampleTable.txt"),
                   "| while read fn ff grp; do echo processing sample ${fn};",
                   fastpPath, 
                   "--qualified_quality_phred", minScore,
                   "-i ${ff} -o", fastpOutputR1,
                   "-w 4 -h ${fn}_fastp.html -j ${fn}_fastp.json -y &&",
                   kallistoPath, "quant -i",
                   file.path(databasePath, database, "index/transcripts.idx"),
                   "-o ${fn}_kallisto -b 100 -t 4",
                   "--single -l", aveFragLen, 
                   "-s", stdFragLen,
                   fastpOutputR1,
                   " && mv ${fn}_kallisto/abundance.tsv ${fn}_kallisto_quant.tsv && ", 
                   " mv ${fn}_kallisto/run_info.json ${fn}_kallisto_log.json && ",
                   "rm -rf ${fn}_kallisto && ",
                   "echo Fastp and Kallisto finish sample ${fn}; done && rm",
                   file.path(userDir, "tmp_clean_*.gz"),
                   "&& echo Fastp and Kallisto have processed all samples;",
                   sep = " ");
    }
    
    strcmd <- paste0("fastp: ", 
                     "\nmininum quality score: ", minScore,
                     "\nkallisto: ", 
                     "\ndatabase: ", list.files(file.path(databasePath, database), pattern = ".gz$"),
                     "\nproject_description: ", des);
    
    str_R <- paste("echo Starting process tables and figures in R $(date);",
                   "Rscript --vanilla",
                   file.path(shellscriptDir, "run_process_kallisto.R"),
                   userDir,
                   "; echo Finish process tables and figures in R $(date);",
                   sep = " ");
    
    str_checMysql1 <- paste("if ! [[ -x ",
                            file.path(userDir, "renameDirFromMysql.sh"),
                            " ]]; then chmod +x ",
                            file.path(userDir, "renameDirFromMysql.sh"),
                            "; fi;",
                            sep = " ");
    
    str_runMysql1 <- paste("bash ",
                           file.path(userDir, "renameDirFromMysql.sh"),
                           ";",
                           sep = " ");
    
    sink(file.path(userDir, "ExecuteRawSeq.sh"));
    cat(str, "\n\n");
    cat(str_R, "\n\n");
    
    cat(str_checMysql1, "\n\n");
    cat(str_runMysql1, "\n\n");
    sink();
    
    sink(file.path(userDir, "analysis_parameters.txt"));
    cat(strcmd, "\n");
    sink();
    
    str_sh <- "#!/bin/sh\n";
    
    str_sh <- paste0(str_sh,
                     "chmod +x ",
                     file.path(userDir, "ExecuteRawSeq.sh; "),
                     "nohup bash ",
                     file.path(userDir, "ExecuteRawSeq.sh > "),
                     file.path(userDir, "slurm_logout.txt 2>&1 &\n"),
                     "echo $! > ",
                     file.path(userDir, "jobid_pid.txt;"));
    
    sink(file.path(userDir, "launcher.sh"));
    cat(str_sh, "\n");
    sink();
    
  }
  
  return(1);
}

SubmitJobDeplexer <- function(userDir, email, database, des, readEnds, shellscriptDir, deplexerExecutable, barcodeRead, barcodeStart, barcodeLength, maxMismatch, trimLeft, sampleIdFile, threads, maxMemory){
  isSlurm <- FALSE
  if(grepl("glassfish", userDir)){
    isSlurm <- TRUE
  }
  
  if(is.null(deplexerExecutable) || deplexerExecutable %in% c("", "NA", "AUTO", "auto")){
    if(file.exists("/home/zgy/NetBeansProjects/")){
      deplexerExecutable <- "/data/glassfish/seq_software/deplexer_source/deplexer"
      isSlurm <- FALSE
    } else if(grepl("glassfish", userDir)){
      deplexerExecutable <- "/data/glassfish/seq_software/deplexer_source/deplexer"
      isSlurm <- TRUE
    } else if(grepl("jeffxia/Dropbox", userDir)) {
      deplexerExecutable <- "/Users/xia/Dropbox/resources/eoa/software/deplexer_source/deplexer"
      isSlurm <- FALSE
    } else {
      deplexerExecutable <- "/home/peng/software/deplexer/deplexer"
      isSlurm <- FALSE
    }
  }
  
  if(readEnds != "pe"){
    stop("DePlexer requires paired-end reads (readEnds == 'pe').")
  }
  
  sampleIdFileNorm <- ifelse(is.null(sampleIdFile), "AUTO", trimws(sampleIdFile))
  if(nchar(sampleIdFileNorm) == 0){
    sampleIdFileNorm <- "AUTO"
  }
  autoSampleId <- tolower(sampleIdFileNorm) == "auto"
  sampleIdPath <- if(grepl("^/", sampleIdFileNorm)) sampleIdFileNorm else file.path(userDir, sampleIdFileNorm)
  sampleIdPathDefault <- file.path(userDir, "Sample_IDs.tsv")
  
  if(isSlurm){
    conf_inf <- paste(
      "#!/bin/bash",
      "#",
      "#SBATCH --job-name=DePlexer_Processing",
      "#",
      "#SBATCH --ntasks=1",
      "#SBATCH --time=14400:00",
      "#SBATCH --mem-per-cpu=5000",
      "#SBATCH --cpus-per-task=4",
      paste0("#SBATCH --mail-user=", email),
      "#SBATCH --mail-type=BEGIN",
      "#SBATCH --mail-type=END",
      "#SBATCH --mail-type=FAIL",
      "#SBATCH --mail-type=REQUEUE",
      "#SBATCH --mail-type=ALL",
      paste0("#SBATCH --output=", file.path(userDir, "slurm_logout.txt")),
      sep = "\n"
    )
    
    if(file.exists("/docker_marker")){
      conf_inf <- paste(
        "#!/bin/bash",
        "#",
        "#SBATCH --job-name=DePlexer_Processing",
        "#",
        "#SBATCH --ntasks=1",
        "#SBATCH --time=28800:00",
        paste0("#SBATCH --output=", file.path(userDir, "slurm_logout.txt")),
        sep = "\n"
      )
    }
  }
  
  sidMode <- if(autoSampleId) "1" else "0"
  sampleIdPathEsc <- gsub("\"", "\\\\\"", sampleIdPath)
  sampleIdPathDefaultEsc <- gsub("\"", "\\\\\"", sampleIdPathDefault)
  sampleIdFileNormEsc <- gsub("\"", "\\\\\"", sampleIdFileNorm)
  userDirEsc <- gsub("\"", "\\\\\"", userDir)
  
  str <- paste(
    sprintf("cd \"%s\" || exit 1;", userDirEsc),
    "cat", file.path(userDir, "sampleTable.txt"),
    "| while read fn ff rf grp; do echo processing sample ${fn};",
    sprintf("if [ -n \"${SLURM_JOB_ID}\" ]; then out_root=\"%s/job_${SLURM_JOB_ID}/output\"; else out_root=\"%s/output\"; fi;", dirname(userDirEsc), userDirEsc),
    "mkdir -p \"${out_root}\";",
    "sample_base=$(basename \"${fn}\");",
    "echo [DePlexer] $(date) sample=${fn} input_R1=${ff} input_R2=${rf};",
    "echo [DePlexer] $(date) sample=${fn} resolving_sample_ids;",
    "sample_dir=$(dirname \"${ff}\");",
    sprintf("user_sid_dir=\"%s\";", userDirEsc),
    "sid_file='';",
    sprintf("if [ \"%s\" = \"1\" ]; then", sidMode),
    "pool=$(echo \"${fn}\" | sed -E 's/^.*\\.([Hh][0-9]+)$/\\1/' | tr '[:upper:]' '[:lower:]');",
    "for d in \"${sample_dir}\" \"${user_sid_dir}\"; do",
    "[ -d \"${d}\" ] || continue;",
    "for p in \"${d}\"/*; do",
    "[ -f \"${p}\" ] || continue;",
    "bn=$(basename \"${p}\");",
    "bn_l=$(printf '%s' \"${bn}\" | tr '[:upper:]' '[:lower:]');",
    "case \"${bn_l}\" in",
    "\"${pool}.sample_ids.tsv\"|\"${pool}.sample_ids.txt\"|\"${pool}.sample_ids.csv\"|\"sample_ids.tsv\"|\"sample_ids.txt\"|\"sample_ids.csv\")",
    "sid_file=\"${p}\"; break;;",
    "esac;",
    "done;",
    "[ -n \"${sid_file}\" ] && break;",
    "done;",
    sprintf("if [ -z \"${sid_file}\" ]; then echo \"ERROR: Missing Sample_IDs file (.tsv/.txt/.csv) for sample ${fn} under ${sample_dir} or %s\"; exit 1; fi;", userDirEsc),
    sprintf("else sid_file=\"%s\";", sampleIdPathEsc),
    "if [ ! -f \"${sid_file}\" ]; then",
    sprintf("sid_file_alt=\"${sample_dir}/%s\";", sampleIdFileNormEsc),
    "if [ -f \"${sid_file_alt}\" ]; then sid_file=\"${sid_file_alt}\"; fi;",
    "if [ ! -f \"${sid_file}\" ]; then",
    "sid_base=$(basename \"${sid_file}\"); sid_base_l=$(printf '%s' \"${sid_base}\" | tr '[:upper:]' '[:lower:]');",
    "for p in \"${sample_dir}\"/*; do [ -f \"${p}\" ] || continue; bn=$(basename \"${p}\"); bn_l=$(printf '%s' \"${bn}\" | tr '[:upper:]' '[:lower:]'); if [ \"${bn_l}\" = \"${sid_base_l}\" ]; then sid_file=\"${p}\"; break; fi; done;",
    "fi;",
    "fi;",
    "if [ ! -f \"${sid_file}\" ]; then echo ERROR: Missing DePlexer sample IDs file ${sid_file}; exit 1; fi;",
    "fi;",
    "echo [DePlexer] $(date) sample=${fn} sample_ids=${sid_file};",
    "start_ts=$(date +%s);",
    "echo [DePlexer] $(date) sample=${fn} deplexer_start;",
    deplexerExecutable,
    "-1 ${ff} -2 ${rf}",
    "-b", barcodeRead,
    "-s", barcodeStart,
    "-l", barcodeLength,
    "-i \"${sid_file}\"",
    "-a", maxMismatch,
    "--trim_left", trimLeft,
    "--prefix ${sample_base}",
    "-o \"${out_root}/${sample_base}_deplexer_out\"",
    "-n", threads,
    "-m", maxMemory,
    "> \"${out_root}/${sample_base}_deplexer_stdout.log\" 2>&1 &",
    "dpid=$!;",
    "echo [DePlexer] $(date) sample=${fn} deplexer_pid=${dpid};",
    "while kill -0 ${dpid} 2>/dev/null; do",
    "if [ -d \"${out_root}/${sample_base}_deplexer_out\" ]; then",
    "n_files=$(find \"${out_root}/${sample_base}_deplexer_out\" -type f | wc -l);",
    "out_sz=$(du -sh \"${out_root}/${sample_base}_deplexer_out\" 2>/dev/null | awk '{print $1}');",
    "else n_files=0; out_sz=0; fi;",
    "echo [DePlexer] $(date) sample=${fn} heartbeat pid=${dpid} out_files=${n_files} out_size=${out_sz};",
    "sleep 30;",
    "done;",
    "wait ${dpid}; status=$?;",
    "end_ts=$(date +%s);",
    "echo [DePlexer] $(date) sample=${fn} deplexer_end status=${status} elapsed_seconds=$((end_ts-start_ts));",
    "if [ ${status} -ne 0 ]; then exit ${status}; fi;",
    "if [ -f \"${out_root}/${sample_base}_deplexer_out/${sample_base}.html\" ]; then",
    "cp -f \"${out_root}/${sample_base}_deplexer_out/${sample_base}.html\" \"${out_root}/${sample_base}_deplexer_report.html\";",
    "echo [DePlexer] $(date) sample=${fn} native_html=${out_root}/${sample_base}_deplexer_out/${sample_base}.html;",
    "echo [DePlexer] $(date) sample=${fn} report_html=${out_root}/${sample_base}_deplexer_report.html;",
    "else echo [DePlexer] $(date) sample=${fn} native_html_missing expected=${out_root}/${sample_base}_deplexer_out/${sample_base}.html; fi;",
    "echo DePlexer finish sample ${fn}; done",
    "&& echo DePlexer have processed all samples",
    "&& echo [DePlexer] $(date) indexing_html_reports;",
    "printf '<html><head><meta charset=\"UTF-8\"><title>DePlexer Reports</title></head><body><h2>DePlexer Reports</h2><ul>' > \"${out_root}/DePlexer_reports_index.html\";",
    "for report in ${out_root}/*_deplexer_report.html; do [ -f \"${report}\" ] && bname=$(basename \"${report}\") && printf '<li><a href=\"%s\">%s</a></li>' \"${bname}\" \"${bname}\" >> \"${out_root}/DePlexer_reports_index.html\"; done;",
    "printf '</ul></body></html>' >> \"${out_root}/DePlexer_reports_index.html\";",
    "echo [DePlexer] $(date) generated_report=${out_root}/DePlexer_reports_index.html;",
    sep = " "
  )
  
  strcmd <- paste0("deplexer: ", deplexerExecutable,
                   "\ndatabase: ", database,
                   "\nbarcode_location: ", barcodeRead,
                   "\nbarcode_start: ", barcodeStart,
                   "\nbarcode_length: ", barcodeLength,
                   "\nbarcode_mismatch: ", maxMismatch,
                   "\ntrim_left: ", trimLeft,
                   "\nsample_ids: ", sampleIdFileNorm,
                   "\nthreads: ", threads,
                   "\nmax_memory: ", maxMemory,
                   "\nproject_description: ", des)
  
  str_checMysql1 <- paste("if ! [[ -x ",
                          file.path(userDir, "renameDirFromMysql.sh"),
                          " ]]; then chmod +x ",
                          file.path(userDir, "renameDirFromMysql.sh"),
                          "; fi;",
                          sep = " ")
  
  str_runMysql1 <- paste("bash ",
                         file.path(userDir, "renameDirFromMysql.sh"),
                         ";",
                         sep = " ")
  
  sink(file.path(userDir, "ExecuteRawSeq.sh"))
  if(isSlurm){
    cat(conf_inf, "\n\n")
  }
  cat(str, "\n\n")
  cat(str_checMysql1, "\n\n")
  cat(str_runMysql1, "\n\n")
  sink()
  
  sink(file.path(userDir, "analysis_parameters.txt"))
  cat(strcmd, "\n")
  sink()
  
  if(!isSlurm){
    str_sh <- "#!/bin/sh\n"
    str_sh <- paste0(str_sh,
                     "chmod +x ",
                     file.path(userDir, "ExecuteRawSeq.sh; "),
                     "nohup bash ",
                     file.path(userDir, "ExecuteRawSeq.sh > "),
                     file.path(userDir, "slurm_logout.txt 2>&1 &\n"),
                     "echo $! > ",
                     file.path(userDir, "jobid_pid.txt;"))
    sink(file.path(userDir, "launcher.sh"))
    cat(str_sh, "\n")
    sink()
  }
  
  return(1)
}

UpdatePcaKls <- function(minReads = 20, projectDirStr, ...) {
  library(data.table)
  library(lattice)
  library(ggplot2)
  library(ggrepel)
  library(vegan)
  library(Cairo)
  library(plyr)
  library(magrittr)
  library(tidyverse)

  mappingReadsDF <- fread(file.path(projectDirStr, "sampleGroup.txt"), header = FALSE)
  names(mappingReadsDF) <-  c("Sample", "Group")

  countfile <- if (file.exists(file.path(projectDirStr, "All_samples_salmon_txi_counts.txt"))) {
    file.path(projectDirStr, "All_samples_salmon_txi_counts.txt")
  } else {
    file.path(projectDirStr, "All_samples_kallisto_txi_counts.txt")
  }
  
  allSamKOAbunDF2 <- fread(countfile) %>% 
    dplyr::slice(2:nrow(.)) %>% 
    column_to_rownames(var = "#NAME")
  allSamKOAbunDF2[] <- lapply(allSamKOAbunDF2, as.numeric)
  allSamKOAbunDF2[allSamKOAbunDF2 < minReads] <- 0
  allSamKOAbunDF2 %<>% filter(rowSums(.) > 0)
  allSamKOAbunDF2 <- t(allSamKOAbunDF2) %>% as.data.frame()
  
  allSamKOAbunDF2 %>% 
    rownames_to_column(var = "sample") %>% 
    fwrite(file.path(projectDirStr, "rareTable.txt"), sep = "\t")
  
  pca <- prcomp(allSamKOAbunDF2, scale. = TRUE)
  pca.res <- as.data.frame(pca$x)
  
  if (length(unique(mappingReadsDF$Group)) > 1) {
    Factor <- mappingReadsDF$Group
  } else {
    Factor <- mappingReadsDF$Sample
  }
  
  pca.res$Group <- Factor
  pca.res$Sample <- rownames(allSamKOAbunDF2)
  
  # Calculate group centroids
  centroids <- aggregate(. ~ Group, data = pca.res[, c("PC1", "PC2", "Group")], mean)
  
  # Merge centroids back to the pca.res dataframe
  pca.res <- merge(pca.res, centroids, by = "Group", suffixes = c("", "_centroid"))
  
  # Calculate the distance to the centroid
  pca.res$distance <- sqrt((pca.res$PC1 - pca.res$PC1_centroid)^2 + (pca.res$PC2 - pca.res$PC2_centroid)^2)
  
  # Identify outliers based on variance threshold (20% here)
  threshold <- 0.2 * mean(pca.res$distance)
  pca.res$outlier <- pca.res$distance > threshold
  
  # Create the PCA plot
  xlim <- range(pca.res$PC1) * 1.1
  ylim <- range(pca.res$PC2) * 1.1
  
  theme.cols <- c("#A7414A","#7A5D1D","#6A8A82", "#CC9B31", "#654242","#A37C27","#281A1A",
                  "#b7c8c3", "#563838", "#DBA3A8", "#435651", "#ceccca", "#282726")
  
  n.factor <- length(unique(Factor))
  fill.cols <- if (n.factor < 14) {
    theme.cols[1:n.factor]
  } else {
    sample(theme.cols, size = n.factor, replace = TRUE)
  }
  
  p <- ggplot(pca.res, aes(x = PC1, y = PC2, color = factor(Group), label = Sample)) +
    geom_point(size = 4) +
    xlim(xlim) + ylim(ylim) +
    labs(color = "Groups") +
    scale_color_manual(values = fill.cols) +
    theme_bw() +
    labs(title = "PCA plot", subtitle = paste0("Min. reads/gene: ", minReads))
  
  if (nrow(pca.res) < 30) {
    p <- p + geom_text_repel(aes(label = Sample), force = 1.5)
  } else {
    p <- p + theme(legend.position = "none")
  }
  
  # Highlight outliers
  p <- p + geom_point(data = pca.res[pca.res$outlier, ], aes(x = PC1, y = PC2), color = "red", shape = 21, size = 4, fill = "red", stroke = 1.5)
  
  # Print the outliers
  outliers <- pca.res[pca.res$outlier, "Sample"]
  cat("Outliers detected:\n")
  
  # Save the plot
  wd <- ifelse(nchar(mappingReadsDF$Sample[1]) < 20, 8, 16)
  ht <- ifelse(nchar(mappingReadsDF$Sample[1]) < 20, 6, 12)
  
  Cairo(file = file.path(projectDirStr, "PCA_sample_similarity.png"), width = wd, height = ht, type = "png", bg = "white", unit = "in", dpi = 300)
  print(p)
  dev.off()
}



UpdateRareKls <- function(step=10, minReads = 20, projectDirStr,
                      facet = NULL, label = NULL, color = NULL, plot = TRUE, 
                      linetype = NULL, se = FALSE, threads = 4,
                      ...){

  library(data.table);
  library(lattice);
  library(ggplot2);
  library(ggrepel);
  library(vegan);
  library(Cairo);
library(plyr);
library(magrittr);
library(tidyverse);
  
  color = "Group";
  label = "Sample";
  #color <- ifelse(is.null(color), "NULL", color);
  linetype <- ifelse(is.null(linetype), "NULL", linetype);
  
  data.src <- fread(file.path(projectDirStr, "metaTable.txt"),
                    header = T); 
  
  countfile = "";
  if(file.exists(file.path(projectDirStr, "All_samples_salmon_txi_counts.txt"))){
    countfile = file.path(projectDirStr, "All_samples_salmon_txi_counts.txt");
  } else {
    countfile = file.path(projectDirStr, "All_samples_kallisto_txi_counts.txt")
  }
  allSamKOAbunDF2 <- fread(countfile) %>% 
    dplyr::slice(2:nrow(.)) %>% 
    column_to_rownames(var = "#NAME");
  allSamKOAbunDF2[] <- lapply(allSamKOAbunDF2, as.numeric);
  allSamKOAbunDF2[] <- lapply(allSamKOAbunDF2, as.integer);
  allSamKOAbunDF2[allSamKOAbunDF2 < minReads] <- 0;
  allSamKOAbunDF2 %<>% 
    filter(rowSums(.) > 0);
  x <- t(allSamKOAbunDF2) %>% as.data.frame();
  x %>% 
    rownames_to_column(var = "sample") %>% 
    fwrite(file.path(projectDirStr, "rareTable.txt"),
           sep = "\t");
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  step_new = floor(max(tot)/as.integer(step))
  rarefun <- function(i) {
    #cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step_new)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- vegan::rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  
  f_n <- paste(data.src, step, "rds", sep = ".");
  
  if(threads > 1){
    out <- parallel::mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE, mc.cores = threads);
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  
  df <- do.call(rbind, out);
  
  # Get sample data
  sdf <- data.src;
  data <- merge(df, sdf, by = "Sample")
  labels <- data.frame(x = tot, y = S, Sample = rownames(x))
  labels <- merge(labels, sdf, by = "Sample")
  
  # Add, any custom-supplied plot-mapped variables
  if (length(color) > 1) {
    data$color <- color
    names(data)[names(data) == "color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  
  if (length(label) > 1) {
    labels$label <- label
    names(labels)[names(labels) == "label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  
  if (length(linetype) > 1) {
    data$linetype <- linetype
    names(data)[names(data) == "linetype"] <- deparse(substitute(linetype))
    linetype <- deparse(substitute(linetype))
  }
  data$Size = as.numeric(data$Size);
  data[["Size_norm"]] <-  data$Size * 100 / as.numeric(data$ReadsMappingRate);
  labels$x = as.numeric(labels$x);
  labels[["x_norm"]] <-  labels$x * 100 / as.numeric(labels$ReadsMappingRate);
  
  p <- ggplot2::ggplot(data = data,
                       ggplot2::aes_string(x = "Size_norm",
                                           y = ".S",
                                           group = "Sample",
                                           color = color,
                                           linetype = linetype))
  
  p <- p + ggplot2::labs(x = "N. of reads", y = "N. of genes");
  
  if (!is.null(label)) {
    p <- p + ggplot2::geom_text(data = labels,
                                ggplot2::aes_string(x = "x_norm",
                                                    y = "y",
                                                    label = label,
                                                    color = color),
                                size = 4, hjust = 0) +
      scale_x_continuous(expand = c(0, 0, 0.3, 0))
  }
  
  p <- p + ggplot2::geom_line();
  if (se) { ## add standard error if available
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes_string(ymin = ".S - .se",
                                               ymax = ".S + .se",
                                               color = NULL,
                                               fill = color),
                           alpha = 0.2)
  }
  
  if(!is.null(facet)){
    p <- p + facet_wrap(as.formula(paste("~", facet)))
  } 
  
  if(nchar(data.src$Sample[1]) < 20){
    wd = 8; ht = 6;
  } else {
    wd = 16; ht = 12;
  }
  
  Cairo(file=file.path(projectDirStr, "RarefactionGene.png"), 
        width=wd, height=ht, type="png", bg="white", unit="in", dpi=300);
  print(p);
  dev.off();
}


pcaPlotS2f <- function(allSamKOAbunDF2,
                       mappingReadsDF, 
                       minReads = 20,
                       projectDirStr) {
  
  # Perform PCA
  pca <- prcomp(allSamKOAbunDF2, scale. = TRUE)
  pca.res <- as.data.frame(pca$x)
  
  # Identify groups
  if(length(unique(mappingReadsDF$Group)) > 1){
    Factor <- mappingReadsDF$Group
  } else {
    Factor <- mappingReadsDF$Sample
  }
  
  pca.res$Group <- Factor
  pca.res$Sample <- row.names(allSamKOAbunDF2)
  
  # Calculate group centroids
  centroids <- aggregate(. ~ Group, data = pca.res[, c("PC1", "PC2", "Group")], mean)
  
  # Merge centroids back to the pca.res dataframe
  pca.res <- merge(pca.res, centroids, by = "Group", suffixes = c("", "_centroid"))
  
  # Calculate the distance to the centroid
  pca.res$distance <- sqrt((pca.res$PC1 - pca.res$PC1_centroid)^2 + (pca.res$PC2 - pca.res$PC2_centroid)^2)
  
  # Identify outliers based on variance threshold (20% here)
  threshold <- 0.2 * mean(pca.res$distance)
  pca.res$outlier <- pca.res$distance > threshold
  
  # Create the PCA plot
  xlim <- range(pca.res$PC1) * 1.1
  ylim <- range(pca.res$PC2) * 1.1
  
  theme.cols <- c("#A7414A","#7A5D1D","#6A8A82", "#CC9B31", "#654242","#A37C27","#281A1A",
                  "#b7c8c3", "#563838", "#DBA3A8", "#435651", "#ceccca", "#282726")
  
  n.factor <- length(unique(Factor))
  if (n.factor < 14) {
    fill.cols <- theme.cols[1:n.factor]
  } else {
    fill.cols <- sample(theme.cols, size = n.factor, replace = TRUE)
  }
  
  p <- ggplot(pca.res, aes(x = PC1, y = PC2, color = factor(Group), label = Sample)) +
    geom_point(size = 4) +
    xlim(xlim) + ylim(ylim) +
    labs(color = "Groups") +
    scale_color_manual(values = fill.cols) +
    theme_bw() +
    labs(title = "PCA plot", subtitle = paste0("Min reads per gene/ortholog: ", minReads))
  
  if (nrow(pca.res) < 30) {
    p <- p + geom_text_repel(aes(label = Sample), force = 1.5)
  }
  
  # Highlight outliers
  p <- p + geom_point(data = pca.res[pca.res$outlier, ], aes(x = PC1, y = PC2), color = "red", shape = 21, size = 4, fill = "red", stroke = 1.5)
  
  # Print the outliers
  outliers <- pca.res[pca.res$outlier, "Sample"]
  cat("Outliers detected:\n")
  print(outliers)
  
  return(p)
}


UpdatePcaS2f <- function(minReads = 20,
                          projectDirStr,
                          ...){
  library(data.table);
  library(lattice);
  library(ggplot2);
  library(ggrepel);
  library(vegan);
  library(Cairo);
library(plyr);
library(magrittr);
library(tidyverse);
  
  metaTable <-  fread(file.path(projectDirStr, "metaTable.txt"),
                      header = T,
                      sep = "\t");
  
  rareTable <- fread(file.path(projectDirStr, "rareTable.txt"),
                     header = T,
                     sep = "\t",
                     nThread = 4);
  rareTable <- as.data.frame(rareTable);
  row.names(rareTable) <- rareTable$V1;
  rareTable <- rareTable[, 2:ncol(rareTable)];
  
  rareTableMinReads <- rareTable;
  rareTableMinReads[rareTableMinReads < minReads] <- 0;
  rareTableMinReads <- rareTableMinReads[, colSums(rareTableMinReads) > 0];
  metaTableMinReads <- metaTable[which(metaTable$Sample == row.names(rareTableMinReads)), ];
  
  #for PCA
  Cairo(file=file.path(projectDirStr, "PCA_sample_similarity.png"), width=8, height=6, type="png", bg="white", unit="in", dpi=300);
  p = pcaPlotS2f(rareTableMinReads, metaTableMinReads, minReads, projectDirStr);
  print(p);
  dev.off();
}


ggRareS2f <- function(s2fObj, data.src, 
                      facet = NULL, label = NULL, color = NULL, plot = TRUE, 
                      linetype = NULL, se = FALSE, step=10, threads = 4,
                      minReads = 20,
                      ...){
  
  color <- ifelse(is.null(color), "NULL", color);
  linetype <- ifelse(is.null(linetype), "NULL", linetype);
  
  x <- s2fObj;
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  step_new = floor(max(tot)/as.integer(step))
  rarefun <- function(i) {
    #cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step_new)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- vegan::rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  
  f_n <- paste(data.src, step, "rds", sep = ".");
  
  if(threads > 1){
    out <- parallel::mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE, mc.cores = threads);
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  
  df <- do.call(rbind, out);
  
  # Get sample data
  sdf <- data.src;
  data <- merge(df, sdf, by = "Sample")
  labels <- data.frame(x = tot, y = S, Sample = rownames(x))
  labels <- merge(labels, sdf, by = "Sample")
  
  # Add, any custom-supplied plot-mapped variables
  if (length(color) > 1) {
    data$color <- color
    names(data)[names(data) == "color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  
  if (length(label) > 1) {
    labels$label <- label
    names(labels)[names(labels) == "label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  
  if (length(linetype) > 1) {
    data$linetype <- linetype
    names(data)[names(data) == "linetype"] <- deparse(substitute(linetype))
    linetype <- deparse(substitute(linetype))
  }
  data$Size = as.numeric(data$Size);
  data[["Size_norm"]] <-  data$Size * 100 / as.numeric(data$percentageReads);
  labels$x = as.numeric(labels$x);
  labels[["x_norm"]] <-  labels$x * 100 / as.numeric(labels$percentageReads);
  
  p <- ggplot2::ggplot(data = data,
                       ggplot2::aes_string(x = "Size_norm",
                                           y = ".S",
                                           group = "Sample",
                                           color = color,
                                           linetype = linetype))
  
  p <- p + ggplot2::labs(x = "N. of reads", y = "N. of Orthologs");
  
  if (!is.null(label)) {
    p <- p + ggplot2::geom_text(data = labels,
                                ggplot2::aes_string(x = "x_norm",
                                                    y = "y",
                                                    label = label,
                                                    color = color),
                                size = 4, hjust = 0) +
      scale_x_continuous(expand = c(0, 0, 0.3, 0))
  }
  
  p <- p + ggplot2::geom_line();
  if (se) { ## add standard error if available
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes_string(ymin = ".S - .se",
                                               ymax = ".S + .se",
                                               color = NULL,
                                               fill = color),
                           alpha = 0.2)
  }
  
  if(!is.null(facet)){
    p <- p + facet_wrap(as.formula(paste("~", facet)))
  } 
  
  invisible(p)
}

UpdateRareS2f <- function(step = 30,
                           minReads = 20,
                           percentageCutoffOri = 75,
                           projectDirStr,
                           ...){
  library(data.table);
  library(lattice);
  library(ggplot2);
  library(ggrepel);
  library(vegan);
  library(Cairo);
library(plyr);
library(magrittr);
library(tidyverse);
  
  percentageCutoffOri = as.numeric(percentageCutoffOri) / 100;
  
  metaTable <- fread(file.path(projectDirStr, "metaTable.txt"),
                     header = T,
                     sep = "\t");
  
  rareTable <- fread(file.path(projectDirStr, "rareTable.txt"),
                     header = T,
                     sep = "\t",
                     nThread = 4);
  rareTable <- as.data.frame(rareTable);
  row.names(rareTable) <- rareTable$V1;
  rareTable <- rareTable[, 2:ncol(rareTable)];
  
  rareTableMinReads <- rareTable;
  rareTableMinReads[rareTableMinReads < minReads] <- 0;
  rareTableMinReads <- rareTableMinReads[, colSums(rareTableMinReads) > 0];
  metaTableMinReads <- metaTable[which(metaTable$Sample == row.names(rareTableMinReads)), ];
  
  #for rarefaction curve
  Cairo(file=file.path(projectDirStr, "RarefactionOrtho.png"), 
        width=8, height=6, type="png", bg="white", unit="in", dpi=300);
  p <- ggRareS2f(s2fObj = rareTableMinReads,
                 data.src = metaTableMinReads,
                 color = "Group",
                 label = "Sample",
                 facet = NULL,
                 linetype = NULL,
                 minReads = minReads,
                 se = FALSE,  # this is not to meaningful
                 step = step, threads = 4);
  print(p)
  dev.off();
  
}
