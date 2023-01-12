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
  gginnards,
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
    # bootstrap
    geom_label_repel(
      aes(x = branch, label = bootstrap, fill = bootstrap),
      #fill = 'white',
      color = "black",
      size = 2,
      hjust =  0.4
    ) +
    scale_fill_gradient(low = "lightblue", high = "lightyellow") +
    theme(legend.position = "none") +
    # node number
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
    style = "height: 100vh; overflow-y: auto;",
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
        "Contribute User Tree File" = "file",
        "Search Lineage in the PR2 Database" = "PR2"
      ),
    ),
    helpText(
      "If you do not have a tree file you can click on the second option to search for the sequences of the target lineage in the PR2 database."
    ),
    
    # User Tree File
    conditionalPanel(
      condition = "input.start == 'file'",
      p(style = "border-bottom: 1px dotted; border-color:#0063B1"),
    ),
    
    # PR2 Pipeline
    conditionalPanel(
      condition = "input.start == 'PR2'",
      p(style = "border-bottom: 1px dotted; border-color:#0063B1"),
      helpText(
        "A) Choose your taxonomic category and lineage group, then click 'Generate'
        to create your fasta sequences as 'pr2_CLADE.fa'."
      ),
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
      fluidRow(
        column(
          width = 8,
          textInput(
            inputId = "clade",
            label = "Lineage Group of Interest:",
            value = "Suessiales"
          )
        ),
        column(
          width = 4,
          style = "margin-top: 25px;",
          actionButton(
            inputId = "pr2clade",
            label = "Generate",
            style = "color: #F9FBFC; background-color: #0063B1; border-color: #0063B1",
            icon = icon("file-export")
          )
        )
      ),
      
      helpText("B) Select your Outgroup Fasta File 'outgroup.fa'. (Optional)"),
      fluidRow(column(
        width = 8,
        fileInput(inputId = "outgroup", label = "Choose Outgroup File:")
      )),
      
      helpText(
        "C) Select the newly generated file from B) named 'pr2_CLADE.fa' and click 'Run Pipeline'."
      ),
      fluidRow(
        column(width = 8,
               fileInput(inputId = "seq", label = "Choose FA File:")),
        column(
          width = 4,
          style = "margin-top: 25px;",
          actionButton(
            inputId = "seqbut",
            label = "Run Pipeline",
            style = "color: #F9FBFC; background-color: #0063B1; border-color: #0063B1",
            icon = icon("circle-play")
          )
        )
      ),
      
      helpText(
        "D) Finally, select the file named 'RAxML_bipartitionsBranchLabels.tre'."
      ),
    ),
    
    # input tree file
    fluidRow(
      column(width = 8,
             fileInput(inputId = "tre", label = "Choose Tree File:")),
      column(
        width = 4,
        style = "margin-top: 25px;",
        actionButton(
          inputId = "treUser",
          label = "Plot Tree",
          style = "color: #F9FBFC; background-color: #0063B1; border-color: #0063B1",
          icon = icon("tree")
        )
      )
    ),
    p(style = "border-bottom: 1px solid"),
    
    ## Section 2
    h3("2) Phylogenetic Tree Editing"),
    helpText("Rerooting must be done before other tree manipulations."),
    
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
          label = "Branch Number:",
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
             numericInput("val_f1", "Node Number:", 87)),
      column(width = 4,
             numericInput("val_f2", "Node Number:", 72)),
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
      "The file must contain the old branch name, atab-separation, and the
               new branch name."
    ),
    fluidRow(
      column(width = 8, fileInput(inputId = "refile", label = "Choose File:")),
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
    p(style = "border-bottom: 1px solid"),
    
    ## Section 4
    h3("4a) Remove Taxa"),
    
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
    p(style = "border-bottom: 1px solid"),
    
    h3("4b) Generate Modified Tree"),
    
    # Modify
    div(h4(em("Modify:")), style = "color:#0063B1"),
    helpText(
      "The modified file has been generated in your directory with the name 'pr2_CLADE_modify.fa'. Outgroup file is optional."
    ),
    fluidRow(
      column(
        width = 8,
        fileInput(inputId = "outgroup2", label = "Choose Outgroup File:")
      ),
      column(
        width = 8,
        fileInput(inputId = "prmod", label = "Choose File with Modifications:")
      ),
      column(
        width = 4,
        style = "margin-top: 25px;",
        actionButton(
          inputId = "removed",
          label = "Run Pipeline",
          icon = icon("circle-play"),
          style = "color: #F9FBFC; background-color: #0063B1; border-color: #0063B1"
        )
      )
    ),
    helpText(
      "To further edit the tree, once complete, return to Step 1) and use the 'RAxML_bipartitionsBranchLabels_modified.tre' as your user input tree."
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
  
  # observeEvent(input$start,
  #               ChartOrder(PutChartOnTop("tree", ChartOrder())))
  
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
  
  # Input SelectBox and TextInput
  mydf <- reactive({
    pr2 <- pr2_database()
    
    group <- switch(
      input$tax,
      "Domain"  = pr2 %>% dplyr::filter(domain  == input$clade) %>% dplyr::select(genbank_accession, sequence_length, sequence),
      "Kingdom" = pr2 %>% dplyr::filter(kingdom == input$clade) %>% dplyr::select(genbank_accession, sequence_length, sequence),
      "Phylum"  = pr2 %>% dplyr::filter(phylum  == input$clade) %>% dplyr::select(genbank_accession, sequence_length, sequence),
      "Class"   = pr2 %>% dplyr::filter(class   == input$clade) %>% dplyr::select(genbank_accession, sequence_length, sequence),
      "Order"   = pr2 %>% dplyr::filter(order   == input$clade) %>% dplyr::select(genbank_accession, sequence_length, sequence),
      "Family"  = pr2 %>% dplyr::filter(family  == input$clade) %>% dplyr::select(genbank_accession, sequence_length, sequence),
      "Genus"   = pr2 %>% dplyr::filter(genus   == input$clade) %>% dplyr::select(genbank_accession, sequence_length, sequence),
      "Species" = pr2 %>% dplyr::filter(species == input$clade) %>% dplyr::select(genbank_accession, sequence_length, sequence)
    )
    
    return(group)
  })
  
  # Function: convert p
  taxonomic <- eventReactive(input$pr2clade, {
    shiny::validate(need(input$clade, "Input Taxonomy and Clade Name."))
    show_modal_spinner(spin = "fading-circle",
                       color = "#0063B1",
                       text = "Generating pr2_CLADE.fa in your Home Directory!")
    
    seq_clade(mydf())
    remove_modal_spinner() # remove it when done
  })
  output$pr2 <- renderDataTable({
    taxonomic()
  })
  
  # PR2
  pl <- eventReactive(input$seqbut, {
    shiny::validate(need(input$seq, "Please Select a pr2_CLADE.fa File!"))
    
    show_modal_spinner(spin = "fading-circle",
                       color = "#0063B1",
                       text = "Please wait...")
    # Vsearch
    update_modal_spinner(text = "Running VSEARCH Sort by Length...")
    x <-
      #c(
      paste(
        "vsearch --sortbylength",
        input$seq$datapath,
        "--output CLADE_sort.fa --minseqlength 500 -notrunclabels",
        sep = " "
      )
    system(x)
    
    update_modal_spinner(text = "Running VSEARCH Cluster...")
    system(
      "vsearch --cluster_smallmem CLADE_sort.fa --id 0.97 --centroids CLADE.clustered.fa -uc CLADE.cluster"
    )
    
    # Import FA files to DNAbin objects
    outgroup_file <- input$outgroup$datapath
    if (!is.null(outgroup_file)) {
      clu <- treeio::read.fasta("CLADE.clustered.fa")
      out <- treeio::read.fasta(input$outgroup$datapath)
      
      # Concatenate the files
      file <- insect::join(clu, out)
    }
    else {
      file <- treeio::read.fasta("CLADE.clustered.fa")
    }
    
    # Saving the sequences as a FA file
    cat(file = "CLADE.cluster.fa",
        paste(paste0(">", names(file)),
              sapply(file, paste, collapse =
                       ""),
              sep = "\n"),
        sep = "\n")
    
    # MAFFT
    update_modal_spinner(text = "Running MAFFT Alignment...")
    system("mafft --reorder --auto CLADE.cluster.fa > CLADE_aligned.fa")
    
    # TrimAl
    update_modal_spinner(text = "Running TRIMAL...")
    system("trimal -in CLADE_aligned.fa -out CLADE.trimal.fa -gt 0.3 -st 0.001")
    
    # RAxML
    update_modal_spinner(text = "Running RAxML...")
    system ("rm -f RAxML_*") #raxml complains if previous files are present, so let's clear them
    system(
      "raxmlHPC-PTHREADS-SSE3 -T 4 -m GTRCAT -c 25 -e 0.001 -p 31415 -f a -N 100 -x 02938 -n tre -s CLADE.trimal.fa"
    )
    
    remove_modal_spinner() # remove it when done
    #}
  })
  output$pr2 <- renderPlot({
    pl()
  })
  
  # PRINT TREES
  user_tree <- eventReactive(input$treUser, {
    shiny::validate(need(input$tre, "Please Select a Tree File!"))
    treeR <- treeio::read.raxml(input$tre$datapath)
    return(treeR)
  })
  output$tree <- renderPlot({
    treeplot(user_tree(), "Phylogenetic Tree")
  })
  
  # ROTATE
  rotate_out <- eventReactive(input$rot, {
    tree <- as.treedata(last_plot())
    
    total_nodes <- treeio::Nnode(tree, internal.only = FALSE)
    internal_nodes <- treeio::Nnode(tree, internal.only = FALSE)
    start_node <- total_nodes - internal_nodes
    
    if (isTRUE(input$val_rot > total_nodes)) {
      shiny::validate("Node number out of range!")
    }
    else {
      tree2 <- last_plot()
      
      tree_rot <- ggtree::rotate(tree2, input$val_rot) +
        geom_hilight(node = input$val_rot, fill = "#0466C8") +
        ggtitle("Rotated Phylogenetic Tree")
      
      tree_rot <-
        move_layers(tree_rot, "GeomHilightRect", position = "bottom")
      
      tree_rot
    }
  })
  output$rotate <- renderPlot({
    rotate_out()
  })
  
  # FLIP
  flip_out <- eventReactive(input$flip, {
    tree <- last_plot()
    
    tree_flip <-
      ggtree::flip(tree, input$val_f1, input$val_f2) + geom_hilight(node = input$val_f1, fill = "#0466C8") +
      geom_hilight(node = input$val_f2, fill = "#002855") +
      ggtitle("Flipped Phylogenetic Tree")
    
    tree_flip <-
      move_layers(tree_flip, "GeomHilightRect", position = "bottom")
    
    tree_flip
  })
  output$flip <- renderPlot({
    flip_out()
  })
  
  # REROOT
  root_out <- eventReactive(input$root, {
    tree <- as.treedata(last_plot())
    all_nodes <- treeio::Nnode(tree, internal.only = FALSE)
    
    if (isTRUE(input$val_root > all_nodes)) {
      shiny::validate(paste0(
        "Error: Node number should be less than or equal to ",
        all_nodes,
        "!"
      ))
    }
    else {
      t <- ape::root(tree,
                     node = input$val_root,
                     resolve.root = TRUE)
      tree_root <- treeplot(t, "Rerooted Phylogenetic Tree")
      tree_root
    }
  })
  output$root <- renderPlot({
    root_out()
  })
  
  # REMOVE
  removee <- eventReactive(input$remove, {
    shiny::validate(need(input$delfile, "Input a Remove File!"))
    shiny::validate(need(input$clufile, "Input a Cluster File!"))
    shiny::validate(need(input$prfile, "Input a pr2_CLADE.fa File!"))
    
    show_modal_spinner(spin = "fading-circle",
                       color = "#0063B1",
                       text = "Generating pr2_CLADE_modify.fa File...")
    
    del <- read_excel(input$delfile$datapath, col_names = 'list')
    del <- del$list
    
    cluster <-
      read.delim(
        input$clufile$datapath,
        header = FALSE,
        sep = "\t",
        dec = "."
      )
    
    dataset <- Biostrings::readDNAStringSet(input$prfile$datapath)
    
    new_clu <- cluster[!cluster$V10 %in% del,]
    new_clu <- new_clu[!new_clu$V9 %in% del,]
    new_clu <- new_clu[!duplicated(new_clu$V9),]
    new_clu <- as.vector(new_clu$V9)
    
    # Export PR2 file
    seq <- names(dataset)
    new_data <- dataset[seq %in% new_clu]
    
    # Saving the sequences as a FA file
    Biostrings::writeXStringSet(new_data, "~/pr2_CLADE_modify.fa", width = 80)
    remove_modal_spinner() # remove it when done
  })
  output$pipeline <- renderText({
    removee()
  })
  
  # MODIFY
  pr2mod <- eventReactive(input$removed, {
    shiny::validate(need(input$prmod, "Please Select the pr2_CLADE_modify.fa File!"))
    
    show_modal_spinner(spin = "fading-circle",
                       color = "#0063B1",
                       text = "Please Wait...!!")
    
    # Vsearch
    update_modal_spinner(text = "Running VSEARCH Sort by Length...")
    x <-
      #c(
      paste(
        "vsearch --sortbylength",
        input$prmod$datapath,
        "--output CLADE_modify_sort.fa --minseqlength 500 -notrunclabels",
        sep = " "
      )
    #)
    system(x)
    
    update_modal_spinner(text = "Running VSEARCH Cluster...")
    system(
      "vsearch --cluster_smallmem CLADE_modify_sort.fa --id 0.97 --centroids CLADE_modify.clustered.fa -uc CLADE_modify.cluster"
    )
    
    # Import FA files to DNAbin objects
    outgroup_file <- input$outgroup2$datapath
    if (!is.null(outgroup_file)) {
      clu <- treeio::read.fasta("CLADE.clustered.fa")
      out <- treeio::read.fasta(input$outgroup2$datapath)
      
      # Concatenate the files
      file <- insect::join(clu, out)
    }
    else {
      file <- treeio::read.fasta("CLADE.clustered.fa")
    }
    
    # Saving the sequences as a FA file
    cat(file = "CLADE_modify.cluster.fa",
        paste(paste0(">", names(file)),
              sapply(file, paste, collapse =
                       ""),
              sep = "\n"),
        sep = "\n")
    
    # MAFFT
    update_modal_spinner(text = "Running MAFFT Alignment...")
    system("mafft --reorder --auto CLADE_modify.cluster.fa > CLADE_modify_aligned.fa")
    
    # TrimAl
    update_modal_spinner(text = "Running TRIMAL...")
    system("trimal -in CLADE_modify_aligned.fa -out CLADE_modify.trimal.fa -gt 0.3 -st 0.001")
    
    # RAxML
    update_modal_spinner(text = "Running RAxML...")
    system ("rm -f RAxML_*") #raxml complains if previous files are present, so let's clear them
    system(
      "raxmlHPC-PTHREADS-SSE3 -T 4 -m GTRCAT -c 25 -e 0.001 -p 31415 -f a -N 100 -x 02938 -n tre -s CLADE_modify.trimal.fa"
    )
    system(
      "mv RAxML_bipartitionsBranchLabels.tre RAxML_bipartitionsBranchLabels_modified.tre"
    )
    
    remove_modal_spinner()
    
    # This is essentially the user input tree from Part 1)
    treeR <-
      treeio::read.raxml("RAxML_bipartitionsBranchLabels_modified.tre")
    
    # Visualize tree
    viz <-
      treeplot(treeR, "Phylogenetic Tree with Modified/Removed Branches")
    
    return(viz)
  })
  output$remove <- renderPlot({
    pr2mod()
  })
  
  # RENAME
  renamee <- eventReactive(input$rename, {
    shiny::validate(need(input$refile, "Input a Rename File!"))
    
    new_name <-
      read.csv(
        input$refile$datapath,
        header = TRUE,
        quote = "\"",
        sep = "\t"
      )
    
    ## this section renames the ggtree graphic names
    #
    tree <- last_plot() + ggtitle("Renamed Phylogenetic Tree")
    # remove the old accession name layer
    tree <- delete_layers(tree, "StatTreeLabel")
    
    # match the accessions to the table replacement names
    viz <- tree %<+% new_name + geom_tiplab(
      aes(label = tax),
      align = TRUE,
      size = 3,
      color = '#609ECF',
      linesize = .3,
      fontface = "bold"
    ) + xlim(NA, .3)
    
    return(viz)
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
        viz <- last_plot()
        viz$data[["label"]] <- viz$data[["tax"]]
        ape::write.tree(as.phylo(viz), file)
      }
    }
  )
}

## APP
shinyApp(ui, server)