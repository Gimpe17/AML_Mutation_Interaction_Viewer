#Load libraries
library(TreeMHN)
library(shiny)
library(DiagrammeR)
library(httr)
library(jsonlite)
library(stringr)
library(shinyBS)



#Load data
load("AML_Data/AML_tree_obj.RData")
load("AML_Data/AML_Theta.RData")



#CODE

#List of drugs (from https://www.cancer.gov/about-cancer/treatment/drugs/leukemia )
drugs <- c(
  "Arsenic Trioxide",
  "Azacitidine",
  "Cyclophosphamide",
  "Cytarabine",
  "Daunorubicin Hydrochloride",
  #"Daunorubicin Hydrochloride and Cytarabine Liposome",
  "Dexamethasone",
  "Doxorubicin Hydrochloride",
  "Enasidenib Mesylate",
  "Gemtuzumab Ozogamicin",
  "Gilteritinib Fumarate",
  "Glasdegib Maleate",
  "Idarubicin Hydrochloride",
  "Ivosidenib",
  "Midostaurin",
  "Mitoxantrone Hydrochloride",
  "Olutasidenib",
  "Pemigatinib",
  "Prednisone",
  "Quizartinib",
  "Rituximab",
  "Thioguanine",
  "Venetoclax",
  "Vincristine Sulfate"
)



# Function to fetch drugs that interact with a gene
get_drugs_interacted <- function(genes) {
  url <- "https://dgidb.org/api/graphql"
  drug_list <- list()
  for (i in genes) {
    query <- sprintf('{
      genes(names: ["%s"]) {
        nodes {
          interactions {
            drug {
              name
            }
            interactionScore
            interactionTypes {
              type
            }
            interactionAttributes {
              name
              value
            }
            sources {
              sourceDbName
            }
          }
        }
      }
    }', i)
    
    # Make the request
    response <- POST(
      url, 
      body = list(query = query), 
      encode = "json",
      content_type_json()
    )
    
    # Check if the request was successful
    if (status_code(response) == 200) {
      data <- content(response, "text", encoding = "UTF-8")
      json_data <- fromJSON(data, flatten = TRUE)
      drug_data <- json_data[["data"]][["genes"]][["nodes"]][["interactions"]][[1]]
      if (length(drug_data) != 0) {
        filtered_data <- subset(drug_data, (tolower(drug.name) %in% tolower(drugs)))
      }
      drug_list[[i]] <- filtered_data
    } else {
      cat(sprintf("Request failed with status code %d\n", status_code(response)))
    }
  }
  return (drug_list)
}

#Creates the drug list
drug_list <- get_drugs_interacted(AML$mutations)



# Function to fetch genes interacted by a drug
get_genes_interacted <- function(drug_name) {
  url <- "https://dgidb.org/api/graphql"
  query <- sprintf('
  {drugs(names: ["%s"]) {
      nodes {
        interactions {
          gene {
            name
  }}}}}', drug_name)
  
  # Make the request
  response <- POST(
    url, 
    body = list(query = query), 
    encode = "json",
    content_type_json()
  )
  
  # Check if the request was successful
  if (status_code(response) == 200) {
    data <- content(response, "text", encoding = "UTF-8")
    json_data <- fromJSON(data, flatten = TRUE)
    genes <- json_data[["data"]][["drugs"]][["nodes"]][["interactions"]][[1]][["gene.name"]]
    #Keep only the genes that are in our database
    filtered_genes <- intersect(genes, AML$mutations)
    return(filtered_genes)
  } else {
    cat(sprintf("Request failed with status code %d\n", status_code(response)))
  }
}



#Modified plot_tree_df function
plot_tree_df1 <- function (tree_df, mutations, tree_label = NULL, main_node_color = "paleturquoise3", 
                           next_node_color = "thistle", genes_interacted = NULL) 
{
  graph_dot <- "\n  digraph g {\n  labelloc='t';\n  fontname='Arial';\n  fontsize=28;\n  "
  if (is.null(tree_label)) {
    graph_dot <- paste(graph_dot, "label = 'Tree", tree_df$Tree_ID[1], "';")
  }
  else {
    graph_dot <- paste(graph_dot, "label = '", tree_label, "';")
  }
  node_labels <- c("Root", mutations[tree_df$Mutation_ID[-1]])
  drug_labels <- c("none", lapply(drug_list[tree_df$Mutation_ID[-1]], function(df) sort(tolower(df[["drug.name"]]))))
  parent_interacted <- rep(0, nrow(tree_df))
  
  if ("Existing" %in% colnames(tree_df)) {
    for (i in c(1:nrow(tree_df))) {
      if (node_labels[i] %in% genes_interacted) {
        if (tree_df$Existing[i]) {
          graph_dot <- paste(graph_dot, tree_df$Node_ID[i], 
                             "[label = '", node_labels[i], "', fontname='Arial', style=filled, color=", 
                             "salmon3", ", tooltip=\"", paste("Drugs:", paste(drug_labels[[i]], collapse = ", ")), "\"];")
        }
        else {
          if ("Prob" %in% colnames(tree_df)) {
            graph_dot <- paste(graph_dot, tree_df$Node_ID[i], 
                               "[label = '", node_labels[i], "\n", tree_df$Prob[i], 
                               "%", "', fontname='Arial', style=filled, color=", 
                               "palegreen3", ", tooltip=\"", paste("Drugs:", paste(drug_labels[[i]], collapse = ", ")), "\"];")
          }
          else {
            graph_dot <- paste(graph_dot, tree_df$Node_ID[i], 
                               "[label = '", node_labels[i], "', fontname='Arial', style=filled, color=", 
                               "palegreen3", ", tooltip=\"", paste("Drugs:", paste(drug_labels[[i]], collapse = ", ")), "\"];")
          }
        }
      parent_interacted[i] <- 1
      } else if (parent_interacted[tree_df$Parent_ID[i]] == 1) {
        if (tree_df$Existing[i]) {
          graph_dot <- paste(graph_dot, tree_df$Node_ID[i], 
                             "[label = '", node_labels[i], "', fontname='Arial', style=filled, color=", 
                             "lightsalmon", ", tooltip=\"", paste("Drugs:", paste(drug_labels[[i]], collapse = ", ")), "\"];")
        }
        else {
          if ("Prob" %in% colnames(tree_df)) {
            graph_dot <- paste(graph_dot, tree_df$Node_ID[i], 
                               "[label = '", node_labels[i], "\n", tree_df$Prob[i], 
                               "%", "', fontname='Arial', style=filled, color=", 
                               "palegreen1", ", tooltip=\"", paste("Drugs:", paste(drug_labels[[i]], collapse = ", ")), "\"];")
          }
          else {
            graph_dot <- paste(graph_dot, tree_df$Node_ID[i], 
                               "[label = '", node_labels[i], "', fontname='Arial', style=filled, color=", 
                               "palegreen1", ", tooltip=\"", paste("Drugs:", paste(drug_labels[[i]], collapse = ", ")), "\"];")
          }
        }
      parent_interacted[i] <- 1
      } else {
        if (tree_df$Existing[i]) {
          graph_dot <- paste(graph_dot, tree_df$Node_ID[i], 
                             "[label = '", node_labels[i], "', fontname='Arial', style=filled, color=", 
                             main_node_color, ", tooltip=\"", paste("Drugs:", paste(drug_labels[[i]], collapse = ", ")), "\"];")
        }
        else {
          if ("Prob" %in% colnames(tree_df)) {
            graph_dot <- paste(graph_dot, tree_df$Node_ID[i], 
                               "[label = '", node_labels[i], "\n", tree_df$Prob[i], 
                               "%", "', fontname='Arial', style=filled, color=", 
                               next_node_color, ", tooltip=\"", paste("Drugs:", paste(drug_labels[[i]], collapse = ", ")), "\"];")
          }
          else {
            graph_dot <- paste(graph_dot, tree_df$Node_ID[i], 
                               "[label = '", node_labels[i], "', fontname='Arial', style=filled, color=", 
                               next_node_color, ", tooltip=\"", paste("Drugs:", paste(drug_labels[[i]], collapse = ", ")), "\"];")
          }
        }
      }
    }
  } else {
    for (i in c(1:nrow(tree_df))) {
      if (node_labels[i] %in% genes_interacted) {
        graph_dot <- paste(graph_dot, tree_df$Node_ID[i], 
                           "[label = '", node_labels[i], "', fontname='Arial', style=filled, color=", 
                           "salmon3", ", tooltip=\"", paste("Drugs:", paste(drug_labels[[i]], collapse = ", ")), "\"];")
        parent_interacted[i] <- 1
      } else if (parent_interacted[tree_df$Parent_ID[i]] == 1) {
        graph_dot <- paste(graph_dot, tree_df$Node_ID[i], 
                           "[label = '", node_labels[i], "', fontname='Arial', style=filled, color=", 
                           "lightsalmon", ", tooltip=\"", paste("Drugs:", paste(drug_labels[[i]], collapse = ", ")), "\"];")
        parent_interacted[i] <- 1
      } else {
        graph_dot <- paste(graph_dot, tree_df$Node_ID[i], 
                         "[label = '", node_labels[i], "', fontname='Arial', style=filled, color=", 
                         main_node_color, ", tooltip=\"", paste("Drugs:", paste(drug_labels[[i]], collapse = ", ")), "\"];")
      }
    }
  }
  for (i in c(2:nrow(tree_df))) {
    graph_dot <- paste(graph_dot, tree_df$Parent_ID[i], "->", 
                       tree_df$Node_ID[i], ";")
  }
  graph_dot <- paste(graph_dot, "}")
  g <- grViz(graph_dot)
  return(g)
}


plot_next_mutations1 <- function (n, tree_df, Theta, mutations = NULL, tree_label = NULL, 
                                 top_M = 5, genes_interacted = NULL) 
{
  if (is.null(mutations)) {
    mutations <- as.character(seq(1, n))
  }
  else {
    if (length(mutations) != n) {
      stop("The number of mutations doesn't match matrix dimension. Please check again...")
    }
    else if (length(unique(mutations)) != n) {
      stop("Mutation names must be unique. Please check again...")
    }
  }
  tree_df$Existing <- rep(1, nrow(tree_df))
  tree_df$Prob <- rep(0, nrow(tree_df))
  next_mutations <- c()
  next_lambdas <- c()
  next_parents <- c()
  for (i in c(1:nrow(tree_df))) {
    pathway_i <- get_pathway_tree_df(tree_df, i)
    siblings <- setdiff(tree_df$Mutation_ID[tree_df$Parent_ID == 
                                              tree_df$Node_ID[i]], c(0))
    for (j in setdiff(c(1:n), c(pathway_i, siblings))) {
      next_mutations <- c(next_mutations, j)
      next_lambdas <- c(next_lambdas, get_lambda(c(pathway_i, 
                                                   j), Theta))
      next_parents <- c(next_parents, tree_df$Node_ID[i])
    }
  }
  cat("Top", top_M, "most probable mutational events that will happen next:\n")
  probs <- next_lambdas/sum(next_lambdas)
  llr <- log(probs * length(probs))
  top_M_idx <- order(probs, decreasing = TRUE)[c(1:top_M)]
  for (i in c(1:top_M)) {
    idx <- top_M_idx[i]
    tree_df <- rbind(tree_df, data.frame(Patient_ID = tree_df$Patient_ID[1], 
                                         Tree_ID = tree_df$Tree_ID[1], Node_ID = nrow(tree_df) + 
                                           1, Mutation_ID = next_mutations[idx], Parent_ID = next_parents[idx], 
                                         Existing = 0, Prob = round(probs[idx] * 100, 3)))
    node <- get_pathway_tree_df(tree_df, nrow(tree_df))
    node <- mutations[node]
    cat("The next most probable node:", paste(c("Root", node), 
                                              collapse = "->"), "\n")
    cat("Probability:", round(probs[idx] * 100, 3), "%\n")
    cat("Log ratio model vs random:", round(llr[idx], 3), 
        "\n")
  }
  g <- plot_tree_df1(tree_df, mutations, tree_label, genes_interacted = genes_interacted)
  return(g)
}

get_pathway_tree_df <- function (tree_df, idx_i) 
{
  pathway <- c()
  repeat {
    mut <- tree_df$Mutation_ID[idx_i]
    if (mut == 0) {
      break
    }
    pathway <- c(mut, pathway)
    idx_i <- match(tree_df$Parent_ID[idx_i], tree_df$Node_ID)
  }
  return(pathway)
}

get_lambda <- function (node, Theta) 
{
  .Call("_TreeMHN_get_lambda", node, Theta)
}






#SHINY APP

# Define UI
ui <- fluidPage(
  titlePanel("AML Mutation Interaction Viewer"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput(inputId = "Patient_ID", label = "Select a patient ID", choices =c("", AML$patients)),
      numericInput(inputId = "n_predictions", label = "Number of predictions", value = 0, min = 0, step = 1),
      selectInput(inputId = "Drug", label = "Select a drug", choices = NULL),
      uiOutput(outputId = "Genes"),
      uiOutput("info"),
      HTML(
        '<div style="margin-bottom: 5px; margin-top: 5px">
          <strong>Legend</strong>
        </div>
        <div style="display: flex; align-items: center; margin-bottom: 12px;">
          <div style="width: 20px; height: 20px; min-width: 20px; min-height: 20px; border-radius: 50%;background-color: #96cdcd; margin-right: 10px;"></div>
          <div>Mutation not interacting with the drug</div>
        </div>
        <div style="display: flex; align-items: center; margin-bottom: 3px;">
          <div style="width: 20px; height: 20px; min-width: 20px; min-height: 20px; border-radius: 50%; background-color: #cd7054; margin-right: 10px;"></div>
          <div>Mutation directly interacting with the drug</div>
        </div>
        <div style="display: flex; align-items: center; margin-bottom: 3px;">
          <div style="width: 20px; height: 20px; min-width: 20px; min-height: 20px; border-radius: 50%; background-color: #ffa07a; margin-right: 10px;"></div>
          <div>Mutation indirectly interacting with the drug (through parent mutation)</div>
        </div>'
      ),
      conditionalPanel(condition = "input.n_predictions != 0", 
             HTML('<div style="display: flex; align-items: center; margin-bottom: 12px;">
                    <div style="width: 20px; height: 20px; min-width: 20px; min-height: 20px; border-radius: 50%; background-color: #d8bfd8; margin-right: 10px;"></div>
                    <div>Predicted mutation not interacting with the drug</div>
                    </div>
                    <div style="display: flex; align-items: center; margin-bottom: 3px;">
                    <div style="width: 20px; height: 20px; min-width: 20px; min-height: 20px; border-radius: 50%; background-color: #7ccd7c; margin-right: 10px;"></div>
                    <div>Predicted mutation directly interacting with the drug</div>
                    </div>
                    <div style="display: flex; align-items: center; margin-bottom: 3px;">
                    <div style="width: 20px; height: 20px; min-width: 20px; min-height: 20px; border-radius: 50%; background-color: #9aff9a; margin-right: 10px;"></div>
                    <div>Predicted mutation indirectly interacting with the drug (through parent mutation)</div>
                    </div>')          
      ),
      width = 3
    ),
    
    mainPanel(
      uiOutput("no_tree"),
      grVizOutput("tree")
    )
  ),
  
  bsModal("modalExample", tags$b("Interaction attributes"), trigger = "ShowAttributes", size = "large",
          uiOutput("Attributes")
  )
)



# Define server logic
server <- function(input, output, session) {
  
  #Update drug selection
  Drugs <- reactive({
    tree_df <- AML$tree_df[AML$tree_df$Tree_ID == match(input$Patient_ID, AML$tree_labels),]
    return(str_to_title(unique(unlist(lapply(drug_list[tree_df$Mutation_ID[-1]], function(df) sort(tolower(df[["drug.name"]])))))))
  })
  observe({
    updateSelectInput(session, "Drug", choices = c("", Drugs()))
  })
  
  #Output info on mutation clicked
  txt <- reactive({
    req(input$tree_click)
    nodeval <- input$tree_click$nodeValues[[1]]
    return(nodeval)
  })
  
  output$info <- renderUI({
    mutation <- trimws(txt())
    if (mutation != "Root") {
      if (tolower(input$Drug) %in% tolower(drug_list[[mutation]][["drug.name"]])) {
        index <- which(tolower(drug_list[[mutation]][["drug.name"]]) == tolower(input$Drug))
        mutation_info <- paste(
          '<p style="margin-top:5px;"><b>Information</b><br>', 
          "Mutation: ", txt(),
          "<br>Drug selected: ", input$Drug,
          "<br>Interaction type: ", drug_list[[mutation]][["interactionTypes"]][[index]],
          "<br>Interaction score: ", format(round(drug_list[[mutation]][["interactionScore"]][[index]], 5), nsmall = 2),
          #"<br>Interaction attributes:" , drug_list[[mutation]][["interactionAttributes"]][[index]],
          "<br>Sources: ", paste(unlist(drug_list[[mutation]][["sources"]][[index]]), collapse = ", "),
          "<br>", actionButton("ShowAttributes", "Show interaction attributes"),
          "</p>"
        )
      } else {
        mutation_info <- paste(
          '<p style="margin-top:5px;"><b>Information</b><br>',
          "Mutation: ", mutation,
          "<br>Interacting drugs:", paste(sort(tolower(drug_list[[which(AML$mutations == mutation)]][["drug.name"]])), collapse = ", "),
          "<p>"
        )
      }
      HTML(mutation_info)
    }  
  })
  
  observeEvent(input$ShowAttributes, {
    toggleModal(session, "modalExample", "open")
  })
  
  output$Attributes <- renderUI({
    mutation <- trimws(txt())
    index <- which(tolower(drug_list[[mutation]][["drug.name"]]) == tolower(input$Drug))
    if (length(drug_list[[mutation]][["interactionAttributes"]][[index]]) == 0 || is.null(drug_list[[mutation]][["interactionAttributes"]][[index]])) {
      print("No attributes")
    } else {
      tableOutput("Attributes_list")
    }},
  )
  
  output$Attributes_list <- renderTable({
    mutation <- trimws(txt())
    index <- which(tolower(drug_list[[mutation]][["drug.name"]]) == tolower(input$Drug))
    drug_list[[mutation]][["interactionAttributes"]][[index]]},
    colnames = FALSE,
    striped = TRUE
  )
  
  
  #Output Graph
  output$no_tree <- renderUI({
    if (is.null(input$Patient_ID) || input$Patient_ID == "") {
      print(h3("Select a patient"))
    }
  })
  output$tree <- renderGrViz({
    if (input$Patient_ID != "") {
      if (input$Drug == "") {
        genes_interacted <- ""
      } else {
        genes_interacted <- get_genes_interacted(input$Drug)
        output$Genes <- renderUI({HTML(paste("<b>Genes interacted by ", input$Drug, ": </b>", paste(genes_interacted, collapse = ", "), sep = ""))})
      }
      if (input$n_predictions == 0) {
        tree <- plot_tree_df1(
          AML$tree_df[AML$tree_df$Tree_ID == match(input$Patient_ID, AML$tree_labels),],
          AML$mutations,
          tree_label = input$Patient_ID,
          genes_interacted = genes_interacted
        )
      } else {
        tree <- plot_next_mutations1(
          AML$n, 
          AML$tree_df[AML$tree_df$Tree_ID == match(input$Patient_ID, AML$tree_labels),],
          AML_Theta,
          AML$mutations,
          tree_label = input$Patient_ID,
          input$n_predictions,
          genes_interacted = genes_interacted
        )
      }
    }
  })
}


# Run the app
shinyApp(ui = ui, server = server)
