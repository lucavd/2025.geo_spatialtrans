#!/usr/bin/env Rscript

# Script per generare un'immagine sintetica di tessuto per simulazioni di trascrittomica spaziale
# L'immagine generata avrà:
# - Un contorno esterno del tessuto con bordo irregolare
# - Diverse regioni interne rappresentanti strutture tissutali distinte
# - Strutture circolari rappresentanti follicoli/ghiandole
# - Strutture lineari/ramificate rappresentanti vasi
# - Gradienti e pattern di densità

suppressPackageStartupMessages(library(imager))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(sp))
set.seed(42)  # Per riproducibilità

generate_synthetic_tissue <- function(
  output_file = "images/synthetic_tissue.png",
  width = 500,
  height = 500,
  n_regions = 5,
  n_structures = 8,
  n_vessels = 3,
  background_value = 0.95,  # Valore di sfondo (bianco = 1, nero = 0)
  tissue_contrast = 0.6,    # Contrasto generale del tessuto
  smoothing_sigma = 3,      # Livello di smoothing per le strutture
  noise_level = 0.05,       # Livello di rumore casuale
  region_overlap = 0.3      # Quanto le regioni possono sovrapporsi
) {
  # Crea un'immagine vuota
  img <- imager::as.cimg(array(background_value, dim = c(width, height, 1, 1)))
  
  # Funzione per creare una regione circolare con bordo sfumato
  create_blob <- function(x_center, y_center, radius, intensity, sigma = smoothing_sigma) {
    # Crea una griglia per calcolare i valori dei pixel
    x_grid <- matrix(rep(1:width, height), ncol = width, byrow = TRUE)
    y_grid <- matrix(rep(1:height, width), ncol = width)
    
    # Calcola la distanza dal centro
    dist_from_center <- sqrt((x_grid - x_center)^2 + (y_grid - y_center)^2)
    
    # Crea un blob con bordo sfumato
    blob <- intensity * exp(-dist_from_center^2 / (2 * radius^2))
    
    # Applica gaussian blur
    if (sigma > 0) {
      blob_img <- as.cimg(blob)
      blob_img <- isoblur(blob_img, sigma)
      blob <- as.array(blob_img)[,,1,1]
    }
    
    return(blob)
  }
  
  # Funzione per creare una forma irregolare (contorno del tessuto)
  create_irregular_shape <- function(x_center, y_center, base_radius, n_points = 20, irregularity = 0.3) {
    # Genera punti attorno al centro con variazione casuale del raggio
    theta <- seq(0, 2*pi, length.out = n_points)
    radius_variation <- runif(n_points, 1 - irregularity, 1 + irregularity)
    radius <- base_radius * radius_variation
    
    # Converti in coordinate cartesiane
    x_points <- x_center + radius * cos(theta)
    y_points <- y_center + radius * sin(theta)
    
    # Chiudi il poligono
    x_points <- c(x_points, x_points[1])
    y_points <- c(y_points, y_points[1])
    
    return(list(x = x_points, y = y_points))
  }
  
  # Funzione per aggiungere rumore a un'immagine
  add_noise <- function(img, noise_level) {
    noise <- noise_level * (rnorm(width * height) %>% 
                             matrix(nrow = height, ncol = width) %>% 
                             as.cimg())
    return(pmax(pmin(img + noise, 1), 0))  # Limita i valori tra 0 e 1
  }
  
  # Funzione per creare un vaso (struttura lineare curva)
  create_vessel <- function(x_start, y_start, x_end, y_end, width, intensity, n_segments = 10) {
    # Crea un percorso curvo tra i punti di inizio e fine
    dx <- (x_end - x_start) / n_segments
    dy <- (y_end - y_start) / n_segments
    
    # Aggiungi casualità al percorso
    x_points <- numeric(n_segments + 1)
    y_points <- numeric(n_segments + 1)
    
    x_points[1] <- x_start
    y_points[1] <- y_start
    
    # Genera punti intermedi con variazione casuale
    for (i in 2:(n_segments+1)) {
      x_points[i] <- x_points[i-1] + dx + rnorm(1, 0, abs(dx) * 0.2)
      y_points[i] <- y_points[i-1] + dy + rnorm(1, 0, abs(dy) * 0.2)
    }
    
    # Crea un'immagine vuota per il vaso
    vessel_img <- as.cimg(array(0, dim = c(width, height, 1, 1)))
    
    # Disegna segmenti connessi con spessore variabile
    for (i in 1:n_segments) {
      segment_width <- width * (0.7 + 0.3 * (n_segments - i) / n_segments)  # Assottigliamento
      blob <- create_blob(x_points[i], y_points[i], segment_width, intensity, sigma = 1)
      vessel_img <- pmax(vessel_img, as.cimg(blob))
    }
    
    return(vessel_img)
  }
  
  # 1. Crea il contorno del tessuto (forma irregolare che copre gran parte dell'immagine)
  tissue_center_x <- width / 2 + rnorm(1, 0, width / 10)
  tissue_center_y <- height / 2 + rnorm(1, 0, height / 10)
  tissue_radius <- min(width, height) * 0.45  # Copre circa il 90% dell'immagine
  
  # Crea forma irregolare per il tessuto
  tissue_shape <- create_irregular_shape(tissue_center_x, tissue_center_y, tissue_radius, 
                                        n_points = 40, irregularity = 0.2)
  
  # Crea un'immagine con il tessuto principale
  tissue_mask <- imager::as.cimg(array(0, dim = c(width, height, 1, 1)))
  
  # Converti i punti del poligono in una maschera
  polygon_points <- data.frame(x = tissue_shape$x, y = tissue_shape$y)
  for (i in 1:width) {
    for (j in 1:height) {
      if (sp::point.in.polygon(i, j, polygon_points$x, polygon_points$y) > 0) {
        tissue_mask[i, j, 1, 1] <- 1
      }
    }
  }
  
  # Applica blur alla maschera del tessuto
  tissue_mask <- isoblur(tissue_mask, sigma = smoothing_sigma * 2)
  
  # Scala la maschera per avere valori tra 0 e tissue_contrast
  tissue_mask <- tissue_mask * tissue_contrast
  
  # 2. Crea regioni all'interno del tessuto
  # Genera centri casuali delle regioni all'interno del contorno
  region_centers <- list()
  region_intensity <- numeric(n_regions)
  region_radius <- numeric(n_regions)
  
  for (i in 1:n_regions) {
    # Trova un punto casuale all'interno del contorno del tessuto
    found_point <- FALSE
    max_tries <- 100
    tries <- 0
    
    while (!found_point && tries < max_tries) {
      # Genera una posizione casuale vicino al centro del tessuto
      candidate_x <- tissue_center_x + rnorm(1, 0, tissue_radius * 0.6)
      candidate_y <- tissue_center_y + rnorm(1, 0, tissue_radius * 0.6)
      
      # Verifica che sia all'interno del contorno
      if (sp::point.in.polygon(candidate_x, candidate_y, polygon_points$x, polygon_points$y) > 0) {
        found_point <- TRUE
        region_centers[[i]] <- c(candidate_x, candidate_y)
        
        # Assegna un'intensità e un raggio casuali
        region_intensity[i] <- runif(1, 0.3, 0.9) * tissue_contrast
        region_radius[i] <- runif(1, tissue_radius * 0.15, tissue_radius * 0.3)
      }
      
      tries <- tries + 1
    }
  }
  
  # Aggiungi le regioni all'immagine
  regions_img <- imager::as.cimg(array(0, dim = c(width, height, 1, 1)))
  
  for (i in 1:n_regions) {
    if (length(region_centers[[i]]) > 0) {
      blob <- create_blob(region_centers[[i]][1], region_centers[[i]][2], 
                        region_radius[i], region_intensity[i])
      
      # Sovrapponi le regioni
      regions_img <- pmax(regions_img, as.cimg(blob) * (1 - region_overlap)) + 
                     as.cimg(blob) * region_overlap * as.numeric(regions_img > 0)
    }
  }
  
  # 3. Aggiungi strutture circolari (ghiandole/follicoli)
  structures_img <- imager::as.cimg(array(0, dim = c(width, height, 1, 1)))
  
  for (i in 1:n_structures) {
    # Trova un punto casuale all'interno del contorno del tessuto
    found_point <- FALSE
    max_tries <- 100
    tries <- 0
    
    while (!found_point && tries < max_tries) {
      # Genera una posizione casuale vicino al centro del tessuto
      candidate_x <- tissue_center_x + rnorm(1, 0, tissue_radius * 0.7)
      candidate_y <- tissue_center_y + rnorm(1, 0, tissue_radius * 0.7)
      
      # Verifica che sia all'interno del contorno
      if (sp::point.in.polygon(candidate_x, candidate_y, polygon_points$x, polygon_points$y) > 0) {
        found_point <- TRUE
        
        # Crea una struttura circolare
        structure_intensity <- runif(1, 0.4, 0.8) * tissue_contrast
        structure_radius <- runif(1, 10, 30)  # Dimensione variabile
        
        # Aggiungi la struttura
        blob <- create_blob(candidate_x, candidate_y, structure_radius, structure_intensity)
        structures_img <- pmax(structures_img, as.cimg(blob))
      }
      
      tries <- tries + 1
    }
  }
  
  # 4. Aggiungi vasi (strutture lineari)
  vessels_img <- imager::as.cimg(array(0, dim = c(width, height, 1, 1)))
  
  for (i in 1:n_vessels) {
    # Genera punto di inizio e fine casuale
    start_angle <- runif(1, 0, 2*pi)
    end_angle <- start_angle + runif(1, pi/2, 3*pi/2)
    
    # Punto di inizio verso il bordo
    start_x <- tissue_center_x + cos(start_angle) * tissue_radius * 0.8
    start_y <- tissue_center_y + sin(start_angle) * tissue_radius * 0.8
    
    # Punto finale verso il centro
    end_x <- tissue_center_x + cos(end_angle) * tissue_radius * 0.3
    end_y <- tissue_center_y + sin(end_angle) * tissue_radius * 0.3
    
    # Assicurati che i punti siano all'interno dell'immagine
    start_x <- max(1, min(width, start_x))
    start_y <- max(1, min(height, start_y))
    end_x <- max(1, min(width, end_x))
    end_y <- max(1, min(height, end_y))
    
    # Crea il vaso
    vessel_width <- runif(1, 5, 15)
    vessel_intensity <- runif(1, 0.3, 0.7) * tissue_contrast
    vessel <- create_vessel(start_x, start_y, end_x, end_y, vessel_width, vessel_intensity)
    
    # Aggiungi il vaso all'immagine
    vessels_img <- pmax(vessels_img, vessel)
  }
  
  # 5. Combina tutti gli elementi
  # Combina in base allo "stack" logico di tessuti
  
  # Normalizza immagine delle regioni
  regions_img <- regions_img / max(regions_img)
  regions_img <- regions_img * tissue_contrast
  
  # Normalizza immagine delle strutture
  structures_img <- structures_img / max(structures_img)
  structures_img <- structures_img * tissue_contrast
  
  # Normalizza immagine dei vasi
  vessels_img <- vessels_img / max(vessels_img)
  vessels_img <- vessels_img * tissue_contrast
  
  # Combina le immagini con pesi diversi
  combined_img <- (
    tissue_mask * 0.5 +      # Base del tessuto
    regions_img * 0.3 +      # Regioni
    structures_img * 0.4 +   # Strutture
    vessels_img * 0.6        # Vasi (più visibili)
  )
  
  # Normalizza l'immagine finale
  combined_img <- combined_img / max(combined_img)
  
  # Inverti l'immagine per avere sfondo chiaro e tessuto scuro
  combined_img <- 1 - combined_img * (1 - background_value)
  
  # Aggiungi rumore per simulare variazioni biologiche
  combined_img <- add_noise(combined_img, noise_level)
  
  # Applica un leggero blur finale per smussare le transizioni
  combined_img <- isoblur(combined_img, sigma = 1)
  
  # Salva l'immagine
  dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
  save.image(combined_img, file = output_file)
  
  # Visualizza l'immagine
  plot(combined_img, main = "Immagine sintetica di tessuto")
  
  # Ritorna l'immagine
  return(combined_img)
}

# Esempio di utilizzo
img <- generate_synthetic_tissue(
  output_file = "images/synthetic_tissue.png",
  width = 500,
  height = 500,
  n_regions = 5,        # Numero di regioni principali
  n_structures = 8,     # Numero di strutture circolari (ghiandole/follicoli)
  n_vessels = 5,        # Numero di vasi
  background_value = 0.95,  # Sfondo quasi bianco
  tissue_contrast = 0.7,    # Contrasto moderato
  smoothing_sigma = 3,      # Livello di smoothing
  noise_level = 0.03,       # Rumore limitato
  region_overlap = 0.2      # Leggera sovrapposizione tra regioni
)

# Genera un'altra versione con parametri leggermente diversi
img2 <- generate_synthetic_tissue(
  output_file = "images/synthetic_tissue2.png",
  width = 500,
  height = 500,
  n_regions = 7,        # Più regioni
  n_structures = 12,    # Più strutture
  n_vessels = 7,        # Più vasi
  background_value = 0.95,
  tissue_contrast = 0.8,    # Maggiore contrasto
  smoothing_sigma = 2,      # Meno smoothing (bordi più netti)
  noise_level = 0.02,
  region_overlap = 0.3      # Maggiore sovrapposizione
)