## Download data from FIA DataMart

for(state in states){
    if (file.exists(file.path(raw_data_dir, paste0(state, "_PLOT.csv"))) &
        file.exists(file.path(raw_data_dir, paste0(state, "_TREE.csv"))) &
        file.exists(file.path(raw_data_dir, paste0(state, "_COND.csv"))) & !redownload_FIA_data) {
        print("Files already downloaded!")
        next
    } else {
        download.file(paste0("https://apps.fs.usda.gov/fia/datamart/CSV/", state, "_PLOT.csv"), 
                      file.path(raw_data_dir, paste0(state, "_PLOT.csv")), 
                      quiet = FALSE, 
                      mode = "w",
                      cacheOK = TRUE)
        download.file(paste0("https://apps.fs.usda.gov/fia/datamart/CSV/", state, "_TREE.csv"), 
                      file.path(raw_data_dir, paste0(state, "_TREE.csv")), 
                      quiet = FALSE, 
                      mode = "w",
                      cacheOK = TRUE)
        download.file(paste0("https://apps.fs.usda.gov/fia/datamart/CSV/", state, "_COND.csv"), 
                      file.path(raw_data_dir, paste0(state, "_COND.csv")), 
                      quiet = FALSE, 
                      mode = "w",
                      cacheOK = TRUE)
   }
}

if(.Platform$OS.type == "unix") {
    info <- system("ls -l", intern=TRUE)
    info <- c(paste0("Files downloaded ", date() , " with redownloading set to ", 
                     redownload_FIA_data), info)
    write(info, file = 'VERSIONS', append = TRUE)
}
