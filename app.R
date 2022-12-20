# PACKAGES
# librarian installs/updates/loads CRAN, github & Bioconductor all-in-one!
if (!require("librarian"))
  install.packages("librarian")
librarian::shelf(
  ape,
  base,
  Biostrings,
  devtools,
  dplyr,
  ggplot2,
  ggrepel,
  ggtree,
  insect,
  pr2database / pr2database,
  readxl,
  shiny,
  shinybusy,
  tidytree,
  treeio,
  quiet = TRUE
)

## FUNCTIONS
seq_clade <- function(x) {
  seq_clade <- Biostrings::DNAStringSet(x$sequence)
  names(seq_clade) <- paste(x$genbank_accession, sep = "|")
  Biostrings::writeXStringSet(seq_clade, "~/pr2_CLADE.fa", width = 80)
}

treeplot <- function(tree, x) {
  ggtree(tree, size = 0.2) +
    geom_tiplab(
      align = TRUE,
      size = 3,
      color = '#609ECF',
      linesize = .3,
      fontface = "bold"
    ) +
    hexpand(.05) +
    labs(title = x) +
    geom_label(
      aes(x = branch, label = bootstrap),
      fill = 'white',
      color = "black",
      size = 2,
      hjust =  0.4
    ) +
    geom_label2(
      aes(subset = !isTip, label = node),
      fill = '#126782',
      color = "white",
      size = 2.9,
      hjust = -0.2
    )
}

## UI
ui <- fluidPage(sidebarLayout(
  sidebarPanel(
    style = "height: 90vh; overflow-y: auto;",
    span(titlePanel(title = div(
      img(
        src = "tra.png",
        height = 72,
        width = 72
      ), "PR2 Curation"
    )), style = "text-align: center"),
    
    ## Header Information
    p(
      span(
        "PR2 Curation is an application that allows users interested in a particular group of microbial eukaryotes to retrieve all sequences belonging to that group, place those sequences in a phylogenetic tree, and curate taxonomic and environmental information about the group."
      ),
      style = "font-size: 16px; color: #00509D"
    ),
    p(style = "border-bottom: 1px dotted; border-color:#0063B1"),
    
    ## Help Button
    actionButton(
      inputId = "help",
      label = "Help",
      icon = icon("circle-question"),
      onclick = "window.open('https://pr2-database.org/', '_blank')",
      style = "color: #F9FBFC; background-color: #0063B1; border-color: #0063B1"
    ),

    ## Section 1
    radioButtons(
      inputId = "start",
      label = h3("1) Creation of Phylogenetic Tree"),
      choices = c(
        "Contribute my Tree File" = "file",
        "Search Lineage in the PR2 Database" = "PR2"
      ),
    ),
    helpText(
      "If you do not have a tree file you can click on the second option to search for the sequences of the target lineage in the PR2 database."
    ),
    
    # User Tree File
    conditionalPanel(
      condition = "input.start == 'file'", 
      fileInput(
        inputId = "tre",
        label = "Choose TRE File:")
    ),
    
    # PR2 Pipeline
    conditionalPanel(
      condition = "input.start == 'PR2'",
      selectInput(
        inputId = "tax",
        label = "Taxonomic Category:",
        choices = c(
          "Domain",
          "Kingdom",
          "Phylum",
          "Class",
          "Order",
          "Family",
          "Genus",
          "Species"
        ),
        selected = "Order"
      ),
      textInput(
        inputId = "clade",
        label = "Lineage Group of Interest:",
        value = "Suessiales"
      ),
      # Help
      # helpText("You can run the example to build a phylogenetic tree
      #        for the order Suessiales, unicellular organisms of the
      #        superclass Dinoflagellata."),
      helpText("Select the file named 'pr2_CLADE.fa'."),
      fluidRow(
        column(width = 8, fileInput(inputId = "seq", label = "Choose FA File:")),
        column(
          width = 4,
          style = "margin-top: 25px;",
          actionButton(
            "seqbut",
            "Run Pipeline",
            style = "color: #F9FBFC; background-color: #0063B1; border-color: #0063B1",
            icon = icon("circle-play"),
          )
        )
      ),
      helpText("Select the file named 'RAxML_bipartitionsBranchLabels.tre'."),
      fluidRow(
        column(width = 8, fileInput(inputId = "raxml", label = "Choose TRE File:")),
        column(
          width = 4,
          style = "margin-top: 25px;",
          actionButton(
            "trePR2",
            "Plot Tree",
            style = "color: #F9FBFC; background-color: #0063B1; border-color: #0063B1",
            icon = icon("tree")
          )
        )
      ),
    ),
    p(style = "border-bottom: 1px solid"),
    
    ## Section 2
    h3("2) Phylogenetic Tree Editing"),
    
    # Reroot
    div(h4(em("Reroot:")), style = "color:#0063B1"),
    fluidRow(
      column(
        width = 8,
        numericInput(
          inputId = "val_root",
          label = "Node number:",
          value = 70
        )
      ),
      column(
        width = 4,
        style = "margin-top: 25px;",
        actionButton(
          inputId = "root",
          label = "Reroot",
          icon = icon("arrows-split-up-and-left"),
          style = "color: #F9FBFC; background-color: #0063B1; border-color: #0063B1"
        )
      )
    ),
    p(style = "border-bottom: 1px solid"),
    
    ## Section 3
    h3("3) Phylogenetic Tree Modification"),
    
    # Rotate
    div(h4(em("Rotate Node:")), style = "color:#0063B1"),
    fluidRow(
      column(
        width = 8,
        numericInput(
          inputId = "val_rot",
          label = "Branch number:",
          value = 87
        )
      ),
      column(
        width = 4,
        style = "margin-top: 25px;",
        actionButton(
          inputId = "rot",
          label = "Rotate",
          icon = icon("rotate"),
          style = "color: #F9FBFC; background-color: #0063B1; border-color: #0063B1"
        )
      )
    ),
    p(style = "border-bottom: 1px dotted; border-color:#0063B1"),
    
    # Flip
    div(h4(em("Flip Nodes:")), style = "color:#0063B1"),
    fluidRow(
      column(width = 4,
             numericInput("val_f1", "Node number:", 87)),
      column(width = 4,
             numericInput("val_f2", "Node number:", 72)),
      column(
        width = 4,
        style = "margin-top: 25px;",
        actionButton(
          inputId = "flip",
          label = "Flip",
          icon = icon("shuffle"),
          style = "color: #F9FBFC; background-color: #0063B1; border-color: #0063B1"
        )
      ),
    ),
    p(style = "border-bottom: 1px dotted; border-color:#0063B1"),
    
    # Rename
    div(h4(em("Rename Nodes:")), style = "color:#0063B1"),
    helpText(
      "The file must contain the old branch name and tab-separated the
               new branch name."
    ),
    fluidRow(
      column(width = 8, fileInput(inputId = "refile", label = "Choose file:")),
      column(
        width = 4,
        style = "margin-top: 25px;",
        actionButton(
          inputId = "rename",
          label = "Rename",
          icon = icon("file-pen"),
          style = "color: #F9FBFC; background-color: #0063B1; border-color: #0063B1"
        )
      )
    ),
    p(style = "border-bottom: 1px dotted; border-color:#0063B1"),
    
    # Remove
    div(h4(em("Remove:")), style = "color:#0063B1"),
    # helpText("The file must contain the old branch name and tab-separated the
    #          new branch name."),
    fileInput(inputId = "delfile",
              label = "Choose file with the names to remove:"),
    fileInput(inputId = "clufile",
              label = "Choose CLUSTER file:"),
    fluidRow(
      column(width = 8,
             fileInput(inputId = "prfile",
                       label = "Choose FA file:")),
      column(
        width = 4,
        style = "margin-top: 25px;",
        actionButton(
          inputId = "remove",
          label = "Remove",
          icon = icon("trash"),
          style = "color: #F9FBFC; background-color: #0063B1; border-color: #0063B1"
        )
      )
    ),
    p(style = "border-bottom: 1px dotted; border-color:#0063B1"),
    
    # Modify
    div(h4(em("Modify:")), style = "color:#0063B1"),
    helpText(
      "The modified file has been generated in your directory with the name 'pr2_CLADE_modify.fa'."
    ),
    fluidRow(
      column(
        width = 8,
        fileInput(inputId = "prmod", label = "Choose file with modifications:")
      ),
      column(
        width = 4,
        style = "margin-top: 25px;",
        actionButton(
          inputId = "removed",
          label = "Plot Tree",
          icon = icon("tree"),
          style = "color: #F9FBFC; background-color: #0063B1; border-color: #0063B1"
        )
      )
    ),
    p(style = "border-bottom: 1px solid"),
    
    ## Save
    radioButtons(
      inputId = "save",
      label = h3("Download Reference Tree"),
      choices =  list("pdf",
                      "tre"),
      selected = "pdf"
    ),
    
    # Download
    downloadButton(
      "down",
      "Download",
      icon = icon("cloud-download"),
      style = "color: #F9FBFC; background-color: #0063B1; border-color: #0063B1"
    )
  ),
  
  # Main panel
  mainPanel(
    dataTableOutput("pr2"),
    dataTableOutput("pipeline"),
    uiOutput("ListOfCharts") #width = "100px", height = "100px"
  )
))

## SERVER
server <- function(input, output) {
  # See only three plots
  PutChartOnTop <- function(This, ChartList) {
    c(This, ChartList[ChartList != This])
  }
  
  ChartOrder <- reactiveVal(list("tree", "root", "rotate", "flip"))
  
  output$ListOfCharts <- renderUI({
    Order <- ChartOrder()[1:3]
    
    ui <- lapply(Order, plotOutput)
    class(ui) <- c("shiny.tag.list", "list")
    return(ui)
  })
  
  observeEvent(input$start,
               ChartOrder(PutChartOnTop("tree", ChartOrder())))
  
  observeEvent(input$val_root | input$root,
               ChartOrder(PutChartOnTop("root", ChartOrder())))
  
  observeEvent(input$val_rot | input$rot,
               {
                 if (input$rot) {
                   ChartOrder(PutChartOnTop("rotate", ChartOrder())) # add plot on top
                 } else {
                   ChartOrder(ChartOrder()[ChartOrder() != "rotate"]) # filter out plot 3
                 }
               })
  
  observeEvent(input$val_f1 | input$flip,
               {
                 if (input$flip) {
                   ChartOrder(PutChartOnTop("flip", ChartOrder()))
                 } else {
                   ChartOrder(ChartOrder()[ChartOrder() != "flip"])
                 }
               })
  
  observeEvent(input$rename,
               {
                 if (input$rename) {
                   ChartOrder(PutChartOnTop("rename", ChartOrder()))
                 } else {
                   ChartOrder(ChartOrder()[ChartOrder() != "rename"])
                 }
               })
  
  observeEvent(input$removed,
               {
                 if (input$removed) {
                   ChartOrder(PutChartOnTop("remove", ChartOrder()))
                 } else {
                   ChartOrder(ChartOrder()[ChartOrder() != "remove"])
                 }
               })
  
  observeEvent(input$trePR2,
               {
                 if (input$trePR2) {
                   ChartOrder(PutChartOnTop("optdos", ChartOrder()))
                 } else {
                   ChartOrder(ChartOrder()[ChartOrder() != "optdos"])
                 }
               })
  
  # Input SelectBox and TextInput
  mydf <- reactive({
    group <- switch(
      input$tax,
      "Domain"  = pr2 %>% dplyr::filter(domain == input$clade)  %>% dplyr::select(genbank_accession, sequence_length, sequence),
      "Kingdom" = pr2 %>% dplyr::filter(kingdom == input$clade) %>% dplyr::select(genbank_accession, sequence_length, sequence),
      "Phylum"  = pr2 %>% dplyr::filter(phylum == input$clade)  %>% dplyr::select(genbank_accession, sequence_length, sequence),
      "Class"   = pr2 %>% dplyr::filter(class == input$clade)   %>% dplyr::select(genbank_accession, sequence_length, sequence),
      "Order"   = pr2 %>% dplyr::filter(order == input$clade)   %>% dplyr::select(genbank_accession, sequence_length, sequence),
      "Family"  = pr2 %>% dplyr::filter(family == input$clade)  %>% dplyr::select(genbank_accession, sequence_length, sequence),
      "Genus"   = pr2 %>% dplyr::filter(genus == input$clade)   %>% dplyr::select(genbank_accession, sequence_length, sequence),
      "Species" = pr2 %>% dplyr::filter(species == input$clade) %>% dplyr::select(genbank_accession, sequence_length, sequence)
    )
    
    return(group)
  })
  
  # Function: convert p
  taxonomic <- reactive({
    req(input$start)
    if (input$start == "PR2") {
      shiny::validate(need(input$clade, "Input correct taxonomy and name."))
      show_modal_spinner(spin = "fading-circle",
                         color = "#0063B1",
                         text = "Please wait...")
      seq_clade(mydf())
      remove_modal_spinner() # remove it when done
    }
  })
  output$pr2 <- renderDataTable({
    taxonomic()
  })
  
  # PR2
  pl <- reactive({
    if (input$seqbut) {
      if (interactive())
        show_modal_spinner(spin = "fading-circle",
                           color = "#0063B1",
                           text = "Please wait...")
      # Vsearch
      x <-
        c(
          paste(
            "vsearch --sortbylength",
            input$seq$datapath,
            "--output CLADE_sort.fa --minseqlength 500 -notrunclabels",
            sep = " "
          )
        )
      system(x)
      system(
        "vsearch --cluster_smallmem CLADE_sort.fa --id 0.97 --centroids CLADE.clustered.fa -uc CLADE.cluster"
      )
      
      # Import FA files to DNAbin objects
      clu <- treeio::read.fasta("CLADE.clustered.fa")
      out <- treeio::read.fasta("outgroup.fa")
      
      # Concatenate the files
      file <- insect::join(clu, out)
      
      # Saving the sequences as a FA file
      cat(file = "CLADE.cluster.fa",
          paste(
            paste0(">", names(file)),
            sapply(file, paste, collapse =
                     ""),
            sep = "\n"
          ),
          sep = "\n")
      # MAFFT
      system("mafft --reorder --auto CLADE.cluster.fa > CLADE_aligned.fa")
      
      # TrimAl
      system("trimal -in CLADE_aligned.fa -out CLADE.trimal.fa -gt 0.3 -st 0.001")
      
      # RAxML
      system(
        "raxmlHPC-PTHREADS-SSE3 -T 4 -m GTRCAT -c 25 -e 0.001 -p 31415 -f a -N 100 -x 02938 -n tre -s CLADE.trimal.fa"
      )
      
      remove_modal_spinner() # remove it when done
    }
  })
  
  #PRINT TREE
  fileee <- reactive({
    # Attach file
    #req(input$start)
    
    if (isTRUE(input$start == "file")) {
      shiny::validate(need(input$tre, "Please Select a File!"))
      treeR <- treeio::read.raxml(input$tre$datapath)
      
      return(treeR)
    }
    else if (isTRUE(input$start == "PR2")) {
    #if (input$trePR2) {
      shiny::validate(need(input$raxml, "Please Select a File!"))
      treeR <- treeio::read.raxml(input$raxml$datapath)
      
      return(treeR)
    }
  })
  
  output$tree <- renderPlot({
    treeplot(fileee(), "Phylogenetic Tree")
  })
  
  # ROTATE
  rotate_out <- reactive({
    if (input$rot) {
      tree <- fileee()
      total_nodes <- treeio::Nnode(tree, internal.only = FALSE)
      internal_nodes <- treeio::Nnode(tree, internal.only = FALSE)
      start_node <- total_nodes - internal_nodes
      
      if (isTRUE(input$val_rot > total_nodes)) {
        shiny::validate("Node number out of range!")
      }
      else {
        viz <- treeplot(tree, "Phylogenetic Tree")
        
        tree_rot <- ggtree::rotate(viz, input$val_rot) +
          geom_hilight(node = input$val_rot, fill = "#0466C8") +
          ggtitle("Rotated Phylogenetic Tree")
        tree_rot
      }
    }
  })
  
  output$rotate <- renderPlot({
    rotate_out()
  })
  
  # FLIP
  flip_out <- reactive({
    if (input$flip) {
      tree <- fileee()
      viz <- treeplot(tree, "Phylogenetic Tree")
      
      tree_flip <-
        ggtree::flip(viz, input$val_f1, input$val_f2) + geom_hilight(node = input$val_f1, fill = "#0466C8") +
        geom_hilight(node = input$val_f2, fill = "#002855") +
        ggtitle("Flipped Phylogenetic Tree")
      tree_flip
    }
  })
  
  output$flip <- renderPlot({
    flip_out()
  })
  
  # REROOT
  root_out <- eventReactive(input$root, {
    #if (input$root) {
      tree <- fileee()
      all_nodes <- treeio::Nnode(tree, internal.only = FALSE)
      
      if (isTRUE(input$val_root > all_nodes)) {
        shiny::validate(paste0("Error: Node number should be less than or equal to ", all_nodes, "!"))
      }
      else {
        t <- ape::root(tree,
                       node = input$val_root,
                       resolve.root = TRUE)
        tree_root <- treeplot(t, "Rerooted Phylogenetic Tree")
        tree_root
      }
   #}
  })
  
  output$root <- renderPlot({
    root_out()
  })
  
  # REMOVE
  removee <- reactive({
    if (input$remove) {
      shiny::validate(need(input$delfile, "Input a File!"))
      show_modal_spinner(spin = "fading-circle",
                         color = "#0063B1",
                         text = "Generating File...")
      
      del <- read_excel(input$delfile$datapath, col_names = FALSE)
      del <- as.vector(del[, 1])
      cluster <-
        read.delim(
          input$clufile$datapath,
          header = FALSE,
          sep = "\t",
          dec = "."
        )
      dataset <- Biostrings::readDNAStringSet(input$prfile$datapath)
      new_clu <- cluster[!cluster$V10 %in% del, ]
      new_clu <- new_clu[!new_clu$V9 %in% del, ]
      new_clu <- new_clu[!duplicated(new_clu$V9), ]
      new_clu <- as.vector(new_clu$V9)
      
      # Export PR2 file
      seq <- names(dataset)
      new_data <- dataset[seq %in% new_clu]
      
      # Saving the sequences as a FA file
      Biostrings::writeXStringSet(new_data, "~/pr2_CLADE_modify.fa", width = 80)
      remove_modal_spinner() # remove it when done
    }
  })
  
  output$pipeline <- renderText({
    removee()
  })
  
  # MODIFY
  pr2mod <- reactive({
    if (input$removed) {
      shiny::validate(need(input$prmod, "Please Select a File!"))
      if (interactive())
        show_modal_spinner(spin = "fading-circle",
                           color = "#0063B1",
                           text = "Please Wait...")
      
      # Vsearch
      update_modal_spinner(text = "Running VSEARCH Sort...")
      x <-
        c(
          paste(
            "vsearch --sortbylength",
            input$prmod$datapath,
            "--output CLADE_sort2.fa --minseqlength 500 -notrunclabels",
            sep = " "
          )
        )
      system(x)
      
      update_modal_spinner(text = "Running VSEARCH Cluster...")
      system(
        "vsearch --cluster_smallmem CLADE_sort2.fa --id 0.97 --centroids CLADE.clustered2.fa -uc CLADE2.cluster"
      )
      #update_modal_progress(1, text = "Finished Running VSEARCH Cluster...")
      
      # Import FA files to DNAbin objects
      clu <- treeio::read.fasta("CLADE.clustered2.fa")
      out <- treeio::read.fasta("outgroup.fa")
      
      # Concatenate the files
      file <- insect::join(clu, out)
      
      # Saving the sequences as a FA file
      cat(file = "CLADE.cluster2.fa",
          paste(
            paste0(">", names(file)),
            sapply(file, paste, collapse =
                     ""),
            sep = "\n"
          ),
          sep = "\n")
      
      # MAFFT
      update_modal_spinner(text = "Running MAFFT Alignment...")
      system("mafft --reorder --auto CLADE.cluster2.fa > CLADE_aligned2.fa")
      
      
      # TrimAl
      update_modal_spinner(text = "Running TRIMAL...")
      system("trimal -in CLADE_aligned2.fa -out CLADE.trimal2.fa -gt 0.3 -st 0.001")
      
      
      # RAxML
      update_modal_spinner(text = "Finished Running RAxML...")
      system ("rm -f RAxML.*") #raxml complains if previous files are present, so let's clear them
      system(
        "raxmlHPC-PTHREADS-SSE3 -T 4 -m GTRCAT -c 25 -e 0.001 -p 31415 -f a -N 100 -x 02938 -n tre -s CLADE.trimal2.fa"
      )
      
      treeR <-
        treeio::read.raxml("RAxML_bipartitionsBranchLabels.tre")
      
      # Visualize tree
      viz <-
        treeplot(treeR, "Phylogenetic Tree with Modified/Removed Branches")
      remove_modal_spinner() # remove it when done
      return(viz)
    }
  })
  
  output$remove <- renderPlot({
    pr2mod()
  })
  
  # RENAME
  renamee <- reactive({
    # Attach file
    if (input$rename) {
      shiny::validate(need(input$refile, "Input a File!"))
      new_name <-
        read.csv(
          input$refile$datapath,
          header = TRUE,
          quote = "\"",
          sep = "\t"
        )
      # replot <- treeplot(treeR, "Phylogenetic tree with renamed branches")
      #treeR@phylo$tip.label[match(new_name$V1, treeR@phylo$tip.label)] <-
      #  new_name$V2
      
      tree <- fileee()
      tree_renamed = rename_taxa(tree, new_name, genbank_accession, tax)
      viz <-
        treeplot(tree_renamed, "Phylogenetic Tree with Renamed Branches ")
      
      #viz %<+% new_name +
      #  geom_tiplab(aes(label=tax))
      
      #p_rename <-
      #  treeplot(viz, "Phylogenetic Tree with Renamed Branches")
      return(viz)
    }
  })
  
  output$rename <- renderPlot({
    renamee()
  })
  
  # SAVE
  output$down <- downloadHandler(
    filename =  function() {
      paste("phylo_tree", input$save, sep = ".")
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      if (input$save == "pdf") {
        ggsave(
          file,
          last_plot(),
          height = 13,
          width = 16.5,
          limitsize = FALSE
        )
      } else if (input$save == "tre") {
        ape::write.tree(as.phylo(last_plot()), file)
      }
    }
  )
}

## APP
shinyApp(ui, server)