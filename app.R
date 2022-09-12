# R PACKAGES
install.packages("shiny")
install.packages("shinybusy")
install.packages("devtools")
install.packages("ape")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("readxl")
install.packages("insect")

# BIOCONDUCTOR PACKAGE
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pr2database")
BiocManager::install("treeio")
BiocManager::install("Biostrings")
BiocManager::install("ggtree")

# LIBRARIES
library(shiny)
library(dplyr)
library(devtools)
library(ggtree)
library(ggplot2)
library(ggrepel)
library(shinybusy)
library(pr2database)
library(readxl)
library(Biostrings) # To save fasta files
library(base)
library(treeio) # tree manipulation
library(ape)
library(insect)

# GITHUB PACKAGE
devtools::install_github("pr2database/pr2database")

## FUNCTIONS
seq_clade <- function(x){
  seq_clade <- Biostrings::DNAStringSet(x$sequence)
  names(seq_clade) <- paste(x$genbank_accession, sep="|")
  Biostrings::writeXStringSet(seq_clade, "~/pr2_CLADE.fa", width = 80)
}

treeplot <- function(tree, x){
  ggtree(tree, size = 0.2) + 
    geom_tiplab(align=TRUE, size=3, color='#609ECF', linesize=.3, fontface = "bold") + # "plain", "bold", "italic", "bold.italic"
    hexpand(.05) +
    labs(title = x) +
    geom_label(aes(x=branch, label=bootstrap), fill='white', color = "black", size = 2, hjust =  0.4) +
    geom_label2(aes(subset=!isTip, label=node), fill='#126782', color = "white", size = 2.9, hjust = -0.2)
}

## UI
ui <- fluidPage(
  
  sidebarLayout(
    
    sidebarPanel(
      style = "height: 90vh; overflow-y: auto;", 
      span(titlePanel(title=div(img(src="tra.png", height = 72, width = 72), "PR2 curation")), style = "text-align: center"),
      p(span("PR2 curation is an application that allows users interested in a particular group of microbial eukaryotes to retrieve all sequences belonging to that group, place those sequences in a phylogenetic tree, and curate taxonomic and environmental information about the group."), style = "font-size: 16px; color: #00509D"), #9B9D9F
      br(), 

      actionButton(inputId = "help", label = "Help", icon = icon("circle-question"), onclick ="window.open('https://pr2-database.org/')", 
                   style="color: #F9FBFC; background-color: #0063B1; border-color: #0063B1"),
      
      radioButtons(inputId = "start",
                   label = h3("Creation of phylogenetic tree"), 
                   choices = c("Contribute my tree file" = "file", 
                               "Search lineage in the PR2 database" = "PR2"),
                   selected = character(0)),
      
      helpText("If you do not have a tree file you can click on the second option to search for the sequences of the target lineage in the PR2 database. "),
      br(),
      
      conditionalPanel(
        condition = "input.start == 'file'",
        fileInput(inputId = "tre",
                  label = "Choose TRE file:")
        ),
    
      conditionalPanel(
        condition = "input.start == 'PR2'",
        selectInput(inputId = "tax",
                    label = "Taxonomic category:",
                    choices = c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                    selected = "Order"),
        textInput(inputId = "clade",
                  label = "Lineage group of interest:",
                  value = "Suessiales"),
        # Help
        # helpText("You can run the example to build a phylogenetic tree 
        #        for the order Suessiales, unicellular organisms of the 
        #        superclass Dinoflagellata."),
        br(),
        helpText("Select the file named 'pr2_CLADE.fa'."),
        fluidRow(
          column(width = 8, fileInput(inputId = "seq", label = "Choose FA file:")),
          column(width = 4, style = "margin-top: 25px;", actionButton("seqbut", "Run pipeline", style="color: #F9FBFC; background-color: #0063B1; border-color: #0063B1")
          )), 
        
        br(),
        helpText("Select the file named 'RAxML_bipartitionsBranchLabels.tre'."),
        fluidRow(
          column(width = 8, fileInput(inputId = "raxml", label = "Choose TRE file:")),
          column(width = 4, style = "margin-top: 25px;", actionButton("trePR2", "Plot tree", style="color: #F9FBFC; background-color: #0063B1; border-color: #0063B1")
          )),
        ),
      
      
      
      p("_____________________________"),
      
      h3("Phylogenetic tree modifications"),
      
      fluidRow(
        column(width = 8, 
               numericInput(inputId = "val_rot", label = "Branch number:", value = 87)),
        column(width = 4, style = "margin-top: 25px;",
               actionButton(inputId = "rot", label = "Rotate",style="color: #F9FBFC; background-color: #0063B1; border-color: #0063B1")
        )
      ), 
      
     br(),
 
      fluidRow(
        column(width = 4, 
               numericInput("val_f1", "Node number:", 87)
        ),
        column(width = 4, 
               numericInput("val_f2", "Node number:", 72)
        ),
        column(width = 4, style = "margin-top: 25px;",
               actionButton(inputId = "flip", label = "Flip", style="color: #F9FBFC; background-color: #0063B1; border-color: #0063B1")
        )
      ),
      
      p("_____________________________"),
      
      h3("Phylogenetic tree editions"),
     div(h4(em("Reroot")), style = "color:#0063B1"),
     fluidRow(
       column(width = 8, 
              numericInput(inputId = "val_root", label = "Node number:", value = 70)),
       column(width = 4, style = "margin-top: 25px;",
              actionButton(inputId = "root", label = "Reroot", style="color: #F9FBFC; background-color: #0063B1; border-color: #0063B1")
       )
     ), 
      
      br(),
      div(h4(em("Remove")), style = "color:#0063B1"),
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
      column(width = 4, style = "margin-top: 25px;",
      actionButton(inputId = "remove",
                   label = "Remove", 
                   style="color: #F9FBFC; background-color: #0063B1; border-color: #0063B1")
     )
    ),
      br(),
    
      helpText("The modified file has been generated in your directory with the name 'pr2_CLADE_modify.fa'."),
     fluidRow(
       column(width = 8, fileInput(inputId = "prmod", label = "Choose file with modifications:")),
       column(width = 4, style = "margin-top: 25px;", actionButton(inputId = "removed", label = "Plot tree", style="color: #F9FBFC; background-color: #0063B1; border-color: #0063B1")
       )
     ), 
     
     div(h4(em("Rename")), style = "color:#0063B1"),
      helpText("The file must contain the old branch name and tab-separated the 
               new branch name."),
     fluidRow(
       column(width = 8, fileInput(inputId = "refile", label = "Choose file:")),
       column(width = 4, style = "margin-top: 25px;", actionButton(inputId = "rename", label = "Rename", style="color: #F9FBFC; background-color: #0063B1; border-color: #0063B1")
       )
     ), 
    
      p("_____________________________"),
      radioButtons(inputId = "save",
                  label = h3("Download reference tree"),
                  choices =  list("pdf",
                               "tre"),
                  selected = "Phylogenetic tree"),
      downloadButton("down", "Download",  icon = icon("cloud-download"), style="color: #F9FBFC; background-color: #0063B1; border-color: #0063B1")
      ),
    
    # Main panel
    mainPanel(
      dataTableOutput("pr2"), 
      dataTableOutput("pipeline"),
      uiOutput("ListOfCharts") #width = "100px", height = "100px"
    )
  )
)

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
 
  observeEvent(
    input$start,
    ChartOrder(PutChartOnTop("tree", ChartOrder())))
  
  observeEvent(
    input$val_root | input$root,
    ChartOrder(PutChartOnTop("root", ChartOrder())))
  
  observeEvent(
    input$val_rot | input$rot,
    {
      if (input$rot) {
        ChartOrder(PutChartOnTop("rotate", ChartOrder())) # add plot on top
      } else {
        ChartOrder(ChartOrder()[ChartOrder() != "rotate"]) # filter out plot 3
      }
    })
  
  observeEvent(
    input$val_f1 | input$flip,
    {
      if (input$flip) {
        ChartOrder(PutChartOnTop("flip", ChartOrder()))
      } else {
        ChartOrder(ChartOrder()[ChartOrder() != "flip"]) 
      }
    })
  
  observeEvent(
    input$rename,
    {
      if (input$rename) {
        ChartOrder(PutChartOnTop("rename", ChartOrder())) 
      } else {
        ChartOrder(ChartOrder()[ChartOrder() != "rename"]) 
      }
    })
  
  observeEvent(
    input$removed,
    {
      if (input$removed) {
        ChartOrder(PutChartOnTop("remove", ChartOrder())) 
      } else {
        ChartOrder(ChartOrder()[ChartOrder() != "remove"]) 
      }
    })
  
  observeEvent(
    input$trePR2,
    {
      if (input$trePR2) {
        ChartOrder(PutChartOnTop("optdos", ChartOrder())) 
      } else {
        ChartOrder(ChartOrder()[ChartOrder() != "optdos"]) 
      }
    })
  
  # Input SelectBox and TextInput
  mydf <- reactive({
    group <- switch(input$tax,
                    "Domain" = pr2 %>% dplyr::filter(domain == input$clade) %>% dplyr::select(genbank_accession, sequence_length, sequence),
                    "Kingdom" = pr2 %>% dplyr::filter(kingdom == input$clade) %>% dplyr::select(genbank_accession, sequence_length, sequence),
                    "Phylum" = pr2 %>% dplyr::filter(phylum == input$clade) %>% dplyr::select(genbank_accession, sequence_length, sequence),
                    "Class" = pr2 %>% dplyr::filter(class == input$clade) %>% dplyr::select(genbank_accession, sequence_length, sequence),
                    "Order" = pr2 %>% dplyr::filter(order == input$clade) %>% dplyr::select(genbank_accession, sequence_length, sequence),
                    "Family" = pr2 %>% dplyr::filter(family == input$clade) %>% dplyr::select(genbank_accession, sequence_length, sequence),
                    "Genus" = pr2 %>% dplyr::filter(genus == input$clade) %>% dplyr::select(genbank_accession, sequence_length, sequence),
                    "Species" = pr2 %>% dplyr::filter(species == input$clade) %>% dplyr::select(genbank_accession, sequence_length, sequence))
    
    return(group)
  })
  
  # Function: convert p
 taxonomic <- reactive({
   req(input$start)
    if(input$start == "PR2"){
      shiny::validate(
        need(input$clade, "Input correct taxonomy and name.")
      )
      show_modal_spinner(
        spin = "fading-circle",
        color = "#0063B1",
        text = "Please wait..."
      )
        seq_clade(mydf())
        remove_modal_spinner() # remove it when done
    }
  })
  
  output$pr2 <- renderDataTable({
    taxonomic()
  })

  # PR2
  pl <- reactive({
    if(input$seqbut){
      if (interactive())
        show_modal_spinner(
          spin = "fading-circle",
          color = "#0063B1",
          text = "Please wait..."
        )
      # Vsearch
      x <- c(paste("vsearch --sortbylength", input$seq$datapath,"--output CLADE_sort.fa --minseqlength 500 -notrunclabels", sep=" "))
      system(x)
      system("vsearch --cluster_smallmem CLADE_sort.fa --id 0.97 --centroids CLADE.clustered.fa -uc CLADE.cluster")
      
      # Import FA files to DNAbin objects
      clu <- treeio::read.fasta("CLADE.clustered.fa")
      out <- treeio::read.fasta("outgroup.fa")
      
      # Concatenate the files
      file <- insect::join(clu, out)
      
      # Saving the sequences as a FA file
      cat(file="CLADE.cluster.fa", paste(paste0(">",names(file)),
                                          sapply(file, paste, collapse=""), sep="\n"), sep="\n")
      # MAFFT
      system("mafft --reorder --auto CLADE.cluster.fa > CLADE_aligned.fa")
      
      # TrimAl
      system("trimal -in CLADE_aligned.fa -out CLADE.trimal.fa -gt 0.3 -st 0.001")

      # RAxML
      system("raxmlHPC-PTHREADS-SSE3 -T 4 -m GTRCAT -c 25 -e 0.001 -p 31415 -f a -N 100 -x 02938 -n tre -s CLADE.trimal.fa")
      
      remove_modal_spinner() # remove it when done
    }
  })
  
  #PRINT TREE
  fileee <- reactive({
    # Attach file
    req(input$start)
    if(input$start == "file"){
      shiny::validate(
        need(input$tre, "Input a file!")
      )
      treeR <- treeio::read.raxml(input$tre$datapath)
      
      # Visualize tree
      viz <- treeplot(treeR, "Phylogenetic tree")
      return(viz)
    }
  })
  
  output$tree <- renderPlot({
    fileee()
  }) # width = 900, height = 450
  
  ### PR2 ###
  pl2 <- reactive({
    # # Attach file
    # req(input$start)
    if(input$trePR2){
      shiny::validate(
        need(input$raxml, "Input a file!")
      )
      treeR <- treeio::read.raxml(input$raxml$datapath)

      # Visualize tree
      viz <- treeplot(treeR, "Phylogenetic tree")
      return(viz)
    }
  })

  output$optdos <- renderPlot({
    pl2()
  })
  
  # ROTATE
  rotate_out <- reactive({
    if(input$rot){
      tree_rot <- ggtree::rotate(fileee(), input$val_rot) +
        geom_hilight(node = input$val_rot, fill = "#0466C8") +
        ggtitle("Rotated phylogenetic tree")
      tree_rot
    }
  })
  
  output$rotate <- renderPlot({
    rotate_out()
  })
  
  # FLIP
  flip_out <- reactive({
    if(input$flip){
      tree_flip <- ggtree::flip(fileee(), input$val_f1, input$val_f2) + geom_hilight(node = input$val_f1, fill = "#0466C8") + 
        geom_hilight(node = input$val_f2, fill = "#002855") + 
        ggtitle("Flipped phylogenetic tree")
      tree_flip
    }
  })
  output$flip <- renderPlot({
    flip_out()
  })
  
  # ROOT
  root_out <- reactive({
    if(input$root){
      t <- ape::root(treeR, node = input$val_root, resolve.root = TRUE)
      tree_root <- treeplot(t, "Enrooted phylogenetic tree")
      tree_root
    }
  })
  output$root <- renderPlot({
    root_out()
  })
  
  # REMOVE
  removee <- reactive({
    if(input$remove){
      show_modal_spinner(
        spin = "fading-circle",
        color = "#0063B1",
        text = "Generating file..."
      )
      del <- read_excel(input$delfile$datapath, col_names = FALSE)
      del <- as.vector(del[,1])
      cluster <- read.delim(input$clufile$datapath, header = FALSE,  sep = "\t", dec = ".")
      dataset <- Biostrings::readDNAStringSet(input$prfile$datapath)
      new_clu <- cluster[! cluster$V10 %in% del, ]
      new_clu <- new_clu[! new_clu$V9 %in% del, ]
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
  
    pr2mod <- reactive({
    if(input$removed){
      if (interactive())
        show_modal_spinner(
          spin = "fading-circle",
          color = "#0063B1",
          text = "Please wait..."
        )
          
          # Vsearch
          x <- c(paste("vsearch --sortbylength", input$prmod$datapath,"--output CLADE_sort2.fa --minseqlength 500 -notrunclabels", sep=" "))
          system(x)
          system("vsearch --cluster_smallmem CLADE_sort2.fa --id 0.97 --centroids CLADE.clustered2.fa -uc CLADE2.cluster")
          
          # Import FA files to DNAbin objects
          clu <- treeio::read.fasta("CLADE.clustered2.fa")
          out <- treeio::read.fasta("outgroup.fa")
          
          # Concatenate the files
          file <- insect::join(clu, out)
          
          # Saving the sequences as a FA file
          cat(file="CLADE.cluster2.fa", paste(paste0(">",names(file)),
                                              sapply(file, paste, collapse=""), sep="\n"), sep="\n")
          # MAFFT
          system("mafft --reorder --auto CLADE.cluster2.fa > CLADE_aligned2.fa")
          
          # TrimAl
          system("trimal -in CLADE_aligned2.fa -out CLADE.trimal2.fa -gt 0.3 -st 0.001")

          # RAxML
          system("raxmlHPC-PTHREADS-SSE3 -T 4 -m GTRCAT -c 25 -e 0.001 -p 31415 -f a -N 100 -x 02938 -n tre -s CLADE.trimal2.fa")
          
          treeR <- treeio::read.raxml("RAxML_bipartitionsBranchLabels.tre")
          
          # Visualize tree
          viz <- treeplot(treeR, "Phylogenetic tree with removed branches")
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
    if(input$rename){
      new_name <- read.table(input$refile$datapath, quote="\"", comment.char="")
      # replot <- treeplot(treeR, "Phylogenetic tree with renamed branches")
      treeR@phylo$tip.label[match(new_name$V1, treeR@phylo$tip.label)] <- new_name$V2
      p_rename <- treeplot(treeR, "Phylogenetic tree with renamed branches")
      return(p_rename)
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
      if(input$save == "pdf"){
        ggsave(file, last_plot(), height=13, width=16.5, limitsize = FALSE)
      } else if(input$save == "tre"){
        ape::write.tree(as.phylo(last_plot()), file)
      }
    } 
  )
}


## APP
shinyApp(ui, server)
