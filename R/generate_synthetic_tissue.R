#!/usr/bin/env Rscript

# Script per generare un'immagine sintetica di tessuto per simulazioni di trascrittomica spaziale
# L'immagine generata avrà:
# - Un contorno esterno del tessuto con bordo circolare
# - Un centro altamente denso con piccole strutture (tipo follicolo germinativo)
# - Diverse strutture circolari rappresentanti follicoli/ghiandole con bordi ben definiti
# - Strutture lineari rappresentanti vasi
# - Texture di base per simulare cellule diffuse

suppressPackageStartupMessages(library(imager))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(sp))
set.seed(42)  # Per riproducibilità

generate_synthetic_tissue <- function(
  output_file = "images/synthetic_tissue.png",
  width = 1000,            # Immagine più grande per maggiore dettaglio
  height = 1000,
  n_follicles = 3,         # Numero di strutture circolari ben definite
  n_vessels = 4,           # Numero di vasi 
  n_cell_clusters = 2000,  # Numero di piccoli cluster cellulari per texture
  background_value = 0.98, # Sfondo molto chiaro (quasi bianco)
  follicle_density = 0.8,  # Densità interna dei follicoli
  noise_level = 0.1,       # Livello di rumore di base
  cell_size = 4            # Dimensione media delle cellule
) {
  # Crea un'immagine di base con sfondo bianco
  img <- imager::as.cimg(array(1.0, dim = c(width, height, 1, 1)))
  
  # Coordinate del centro dell'immagine
  center_x <- width / 2
  center_y <- height / 2
  
  # Raggio del tessuto (circa 80% dell'immagine)
  tissue_radius <- min(width, height) * 0.4
  
  # Funzione per creare un follicolo (struttura rotonda con bordo e texture interna)
  create_follicle <- function(x, y, radius, border_width = radius * 0.15) {
    # Crea una maschera temporanea alla dimensione dell'immagine
    mask <- imager::as.cimg(array(0, dim = c(width, height, 1, 1)))
    
    # Prepara grid di coordinate
    grid_x <- matrix(rep(1:width, each = height), nrow = height, ncol = width)
    grid_y <- matrix(rep(1:height, times = width), nrow = height, ncol = width)
    
    # Calcola distanza dal centro del follicolo
    dist_from_center <- sqrt((grid_x - x)^2 + (grid_y - y)^2)
    
    # Crea il bordo esterno più scuro 
    outer_border <- (dist_from_center <= radius) & (dist_from_center > radius - border_width)
    # Area interna del follicolo
    inner_area <- dist_from_center <= (radius - border_width)
    
    # Assegna i valori
    mask[outer_border] <- 0.6  # Bordo più scuro (valore più basso = più scuro)
    mask[inner_area] <- 0.75   # Area interna un po' più chiara
    
    # Aggiungi texture all'interno (piccoli punti più scuri)
    for (i in 1:round(radius^2 * follicle_density * 0.3)) {
      # Posizione casuale interna al follicolo (escluso il bordo)
      r <- runif(1, 0, radius - border_width - 3)
      theta <- runif(1, 0, 2*pi)
      
      dot_x <- x + r * cos(theta)
      dot_y <- y + r * sin(theta)
      dot_radius <- runif(1, 1.5, 3.5)
      
      # Calcola i pixel da modificare
      dot_mask <- sqrt((grid_x - dot_x)^2 + (grid_y - dot_y)^2) <= dot_radius
      
      # Aggiungi il punto scuro
      mask[dot_mask & inner_area] <- mask[dot_mask & inner_area] * 0.4
    }
    
    # Applica un leggero blur per smussare
    mask <- isoblur(mask, sigma = 1.2)
    
    return(mask)
  }
  
  # Funzione per creare un vaso (struttura lineare curva)
  create_vessel <- function(x_start, y_start, x_end, y_end, vessel_width = 10, curvature = 0.2) {
    # Crea una maschera vuota
    mask <- imager::as.cimg(array(0, dim = c(width, height, 1, 1)))
    
    # Prepara grid di coordinate (riutilizzabile)
    if (!exists("grid_x") || !exists("grid_y")) {
      grid_x <<- matrix(rep(1:width, each = height), nrow = height, ncol = width)
      grid_y <<- matrix(rep(1:height, times = width), nrow = height, ncol = width)
    }
    
    # Numero di segmenti per il vaso
    n_segments <- 20
    
    # Calcola la curva del vaso con un controllo di curvatura
    # Usiamo le curve di Bezier per creare un percorso naturale
    t <- seq(0, 1, length.out = n_segments)
    
    # Punto di controllo per la curva (perpendiculare alla linea tra inizio e fine)
    dx <- x_end - x_start
    dy <- y_end - y_start
    dist <- sqrt(dx^2 + dy^2)
    
    # Punto di controllo per la curva di Bezier
    ctrl_x <- (x_start + x_end) / 2 - dy * curvature
    ctrl_y <- (y_start + y_end) / 2 + dx * curvature
    
    # Calcola i punti della curva di Bezier
    x_curve <- (1-t)^2 * x_start + 2*(1-t)*t * ctrl_x + t^2 * x_end
    y_curve <- (1-t)^2 * y_start + 2*(1-t)*t * ctrl_y + t^2 * y_end
    
    # Disegna il vaso lungo la curva
    for (i in 1:n_segments) {
      # Varia leggermente la larghezza per un aspetto naturale
      segment_width <- vessel_width * (0.8 + 0.4 * sin(i / n_segments * pi))
      
      # Seleziona i pixel nel raggio del vaso
      vessel_mask <- sqrt((grid_x - x_curve[i])^2 + (grid_y - y_curve[i])^2) <= segment_width
      
      # Aggiunge il segmento al vaso (valore più chiaro verso il centro)
      mask[vessel_mask] <- 0.7
    }
    
    # Applica un leggero blur per smussare
    mask <- isoblur(mask, sigma = 1.5)
    
    return(mask)
  }
  
  # Funzione per creare piccoli cluster cellulari per texture
  create_cell_clusters <- function(n_clusters, radius_range = c(1.5, 3.5), intensity_range = c(0.7, 0.9)) {
    # Crea una maschera vuota
    mask <- imager::as.cimg(array(0, dim = c(width, height, 1, 1)))
    
    # Prepara grid di coordinate se non esistono già
    if (!exists("grid_x") || !exists("grid_y")) {
      grid_x <<- matrix(rep(1:width, each = height), nrow = height, ncol = width)
      grid_y <<- matrix(rep(1:height, times = width), nrow = height, ncol = width)
    }
    
    # Genera tutti i cluster in un colpo solo
    # Posizioni casuali all'interno del tessuto
    r <- runif(n_clusters, 0, tissue_radius * 0.95)
    theta <- runif(n_clusters, 0, 2*pi)
    
    x_coords <- center_x + r * cos(theta)
    y_coords <- center_y + r * sin(theta)
    
    # Raggi e intensità
    radii <- runif(n_clusters, radius_range[1], radius_range[2])
    intensities <- runif(n_clusters, intensity_range[1], intensity_range[2])
    
    # Aggiungi tutti i cluster
    for (i in 1:n_clusters) {
      cell_mask <- sqrt((grid_x - x_coords[i])^2 + (grid_y - y_coords[i])^2) <= radii[i]
      mask[cell_mask] <- intensities[i]
    }
    
    return(mask)
  }
  
  # 1. Crea la maschera base del tessuto (un cerchio con bordo sfumato)
  tissue_mask <- imager::as.cimg(array(0, dim = c(width, height, 1, 1)))
  
  # Prepara grid di coordinate
  grid_x <- matrix(rep(1:width, each = height), nrow = height, ncol = width)
  grid_y <- matrix(rep(1:height, times = width), nrow = height, ncol = width)
  
  # Distanza dal centro
  dist_from_center <- sqrt((grid_x - center_x)^2 + (grid_y - center_y)^2)
  
  # Crea la maschera del tessuto con bordo sfumato
  tissue_mask[dist_from_center <= tissue_radius] <- 0.9
  
  # Sfuma il bordo
  blur_width <- tissue_radius * 0.05
  border <- (dist_from_center > tissue_radius - blur_width) & (dist_from_center <= tissue_radius)
  fade_factor <- 1 - (dist_from_center[border] - (tissue_radius - blur_width)) / blur_width
  tissue_mask[border] <- 0.9 * fade_factor
  
  # 2. Crea un grande cluster centrale (tipo centro germinativo)
  central_cluster_radius <- tissue_radius * 0.3
  central_cluster <- dist_from_center <= central_cluster_radius
  
  # Area del cluster centrale
  central_area <- imager::as.cimg(array(0, dim = c(width, height, 1, 1)))
  central_area[central_cluster] <- 0.6  # Più scuro
  
  # Aggiungi micro-texture nel cluster centrale (più densa)
  central_texture <- create_cell_clusters(n_clusters = n_cell_clusters * 0.4, 
                                          radius_range = c(1, 2.5), 
                                          intensity_range = c(0.3, 0.5))
  
  # Limita la texture all'area centrale
  central_texture <- central_texture * as.numeric(central_cluster)
  
  # 3. Crea i follicoli (strutture circolari con bordo e texture interna)
  follicle_img <- imager::as.cimg(array(0, dim = c(width, height, 1, 1)))
  
  for (i in 1:n_follicles) {
    # Posizione distribuita come satelliti intorno al centro
    angle <- (i-1) * (2*pi / n_follicles) + runif(1, -pi/6, pi/6)
    dist_from_center <- tissue_radius * runif(1, 0.4, 0.7)
    
    # Converti in coordinate cartesiane
    x <- center_x + dist_from_center * cos(angle)
    y <- center_y + dist_from_center * sin(angle)
    
    # Raggio variabile
    follicle_radius <- tissue_radius * runif(1, 0.1, 0.15)
    
    # Crea e aggiungi il follicolo
    follicle <- create_follicle(x, y, follicle_radius)
    follicle_img <- pmax(follicle_img, follicle)
  }
  
  # 4. Crea i vasi (strutture lineari che attraversano il tessuto)
  vessels_img <- imager::as.cimg(array(0, dim = c(width, height, 1, 1)))
  
  for (i in 1:n_vessels) {
    # Punto di partenza sul bordo
    start_angle <- runif(1, 0, 2*pi)
    start_x <- center_x + tissue_radius * 0.9 * cos(start_angle)
    start_y <- center_y + tissue_radius * 0.9 * sin(start_angle)
    
    # Punto finale verso il centro
    angle_diff <- runif(1, pi/4, 3*pi/4) * sample(c(-1, 1), 1)
    end_angle <- start_angle + angle_diff
    
    end_dist <- tissue_radius * runif(1, 0.1, 0.5)
    end_x <- center_x + end_dist * cos(end_angle)
    end_y <- center_y + end_dist * sin(end_angle)
    
    # Larghezza variabile
    vessel_width <- runif(1, 8, 15)
    
    # Crea e aggiungi il vaso
    vessel <- create_vessel(start_x, start_y, end_x, end_y, vessel_width = vessel_width, 
                           curvature = runif(1, 0.1, 0.3))
    vessels_img <- pmax(vessels_img, vessel)
  }
  
  # 5. Crea texture generale di cellule (meno densa fuori dal centro)
  general_texture <- create_cell_clusters(n_clusters = n_cell_clusters * 0.6, 
                                         radius_range = c(1.5, 3), 
                                         intensity_range = c(0.65, 0.85))
  
  # Limita la texture all'area del tessuto
  general_texture <- general_texture * as.numeric(dist_from_center <= tissue_radius)
  
  # 6. Combina tutte le componenti
  # Inizia con il valore di background
  combined_img <- imager::as.cimg(array(background_value, dim = c(width, height, 1, 1)))
  
  # Inverti le maschere (1 = bianco, 0 = nero, ma vogliamo 0 = bianco, 1 = nero per l'immagine finale)
  # Sottrai solo dove la maschera è > 0 per mantenere il background
  
  # Applica la maschera del tessuto
  tissue_effect <- tissue_mask * 0.1  # Effetto leggero
  combined_img <- combined_img - tissue_effect * as.numeric(tissue_mask > 0)
  
  # Applica la texture generale (effetto leggero)
  texture_effect <- general_texture * 0.25
  combined_img <- combined_img - texture_effect * as.numeric(general_texture > 0)
  
  # Applica il cluster centrale (effetto più forte)
  central_effect <- central_area * 0.3 + central_texture * 0.4
  combined_img <- combined_img - central_effect * as.numeric((central_area + central_texture) > 0)
  
  # Applica i follicoli (con priorità, sovrascrivendo altri elementi)
  follicle_effect <- follicle_img * 0.5
  combined_img <- combined_img - follicle_effect * as.numeric(follicle_img > 0)
  
  # Applica i vasi (con priorità alta)
  vessel_effect <- vessels_img * 0.6
  combined_img <- combined_img - vessel_effect * as.numeric(vessels_img > 0)
  
  # 7. Aggiungi rumore di base per texture più naturale
  noise <- noise_level * matrix(rnorm(width * height), nrow = height, ncol = width)
  noise_img <- as.cimg(noise) 
  
  # Limita il rumore e applica solo all'area del tessuto
  noise_img <- noise_img * as.numeric(dist_from_center <= tissue_radius * 1.05)
  combined_img <- combined_img - abs(noise_img) * 0.15
  
  # 8. Assicurati che i valori rimangano nel range corretto [0,1]
  combined_img <- pmax(pmin(combined_img, 1), 0)
  
  # 9. Applica un leggero blur finale per smussare le transizioni
  combined_img <- isoblur(combined_img, sigma = 0.8)
  
  # 10. Salva l'immagine
  dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
  save.image(combined_img, file = output_file)
  
  # Visualizza l'immagine
  plot(combined_img, main = "Immagine sintetica di tessuto")
  
  # Ritorna l'immagine
  return(combined_img)
}

# Genera una prima versione
img1 <- generate_synthetic_tissue(
  output_file = "images/synthetic_tissue1.png",
  width = 1000,
  height = 1000,
  n_follicles = 3,
  n_vessels = 4,
  n_cell_clusters = 2000,
  cell_size = 4
)

# Genera una seconda versione con parametri leggermente diversi
img2 <- generate_synthetic_tissue(
  output_file = "images/synthetic_tissue2.png",
  width = 1000,
  height = 1000,
  n_follicles = 4,      # Più follicoli
  n_vessels = 5,        # Più vasi
  n_cell_clusters = 3000, # Più cellule
  follicle_density = 0.9  # Follicoli più densi
)