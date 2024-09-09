#Mis funciones 

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

library(pacman)

library(pacman)
p_load(dplyr)

#----------------------------------------------------#
#            COnteo de grupod por columna            #
#----------------------------------------------------#
count_group <- function(column) {
  result <- data.frame(Group = names(table(column)), Count = as.vector(table(column)))
  return(result)
}
    
     #USO: count_group(mi_dato$mi columna) 


#----------------------------------------------------#
#            EXTRAER ORs e INTERVALOS DE CONFIANZA   #
#----------------------------------------------------#
generateORICP <- function(modelo) {
  
  estimates <- summary(modelo)$coefficients[,1]
  p_values <- summary(modelo)$coefficients[,4]
  significance <- ifelse(p_values < 0.001, "***", ifelse(p_values < 0.01, "**", ifelse(p_values < 0.05, "*", ifelse(p_values < 0.1, ".", ""))))
  
  OR <- exp(cbind(coef(modelo), confint(modelo)))
  
  df <- data.frame(Estimate = estimates, OR = OR[,1], "2.5%" = OR[,2], 
                   "97.5%" = OR[,3], "p-value" = p_values, "significancia" = significance)
  
  return(df)
}

    #USO: generateORICP(mi_modelo)

#----------------------------------------------------#
#         EXTRAER ORs e INTERVALOS DE CONFIANZA  V2  #
#----------------------------------------------------#
getORICP <- function(modelo) {
  
  estimates <- summary(modelo)$coefficients[,1]
  p_values <- summary(modelo)$coefficients[,4]
  significance <- ifelse(p_values < 0.001, "***", ifelse(p_values < 0.01, "**", ifelse(p_values < 0.05, "*", ifelse(p_values < 0.1, ".", ""))))
  
  OR <- exp(cbind(coef(modelo), confint(modelo)))
  
  df <- data.frame(Estimate = estimates, OR = OR[,1], "2.5%" = OR[,2], 
                   "97.5%" = OR[,3], "p-value" = p_values, "significancia" = significance)
  
  return(df)
}


#--------------------------------------------------#
#           TRaducir de DNA to amino               #
#--------------------------------------------------#

# No aceptar numeros
#Mensaje: No cadena en multiplo de 3

DNAtoAmino_V1 <- function(cadena) {
  
  # Crear el data frame de codones
  codones_df <- data.frame(
    codon = c(
      "TTT", "TTC", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG",
      "ATT", "ATC", "ATA", "ATG", "GTT", "GTC", "GTA", "GTG",
      "TCT", "TCC", "TCA", "TCG", "CCT", "CCC", "CCA", "CCG",
      "ACT", "ACC", "ACA", "ACG", "GCT", "GCC", "GCA", "GCG",
      "TAT", "TAC", "TAA", "TAG", "CAT", "CAC", "CAA", "CAG",
      "AAT", "AAC", "AAA", "AAG", "GAT", "GAC", "GAA", "GAG",
      "TGT", "TGC", "TGA", "TGG", "CGT", "CGC", "CGA", "CGG", "AGT", "AGC", "AGA", "AGG",
      "GGT", "GGC", "GGA", "GGG"
    ),
    aminoacid = c(
      "F", "F", "L", "L", "L", "L", "L", "L",
      "I", "I", "I", "M", "V", "V", "V", "V",
      "S", "S", "S", "S", "P", "P", "P", "P",
      "T", "T", "T", "T", "A", "A", "A", "A",
      "Y", "Y", "-", "-", "H", "H", "Q", "Q",
      "N", "N", "K", "K", "D", "D", "E", "E",
      "C", "C", "-", "W", "R", "R", "R", "R", "S", "S", "R", "R",
      "G", "G", "G", "G"
    )
  )
  
  # Función interna para traducir una sola cadena
  traducirCadena <- function(cadena) {
    if(nchar(cadena) < 1) {
      return("-") # o return("") si prefieres una cadena vacía en vez de "-"
    }
    
    aminoacidos <- c()
    
    for (i in seq(1, nchar(cadena), 3)) {
      triplete <- substr(cadena, i, i + 2)
      
      # Buscar el triplete en el data frame de codones
      aa <- codones_df$aminoacid[which(codones_df$codon == triplete)]
      
      # Si el triplete no se encontró o no tiene 3 bases, se asigna "-"
      if (length(aa) == 0 || nchar(triplete) < 3) {
        aa <- "-"
      }
      
      aminoacidos <- c(aminoacidos, aa)
    }
    
    return(paste(aminoacidos, collapse = ""))
  }
  
  # Si 'cadena' es un vector, aplicar 'traducirCadena' a cada elemento
  if(is.vector(cadena)) {
    res <- sapply(cadena, traducirCadena)
    names(res) <- NULL  # quitar nombres del vector
    return(res)
  } else {
    return(traducirCadena(cadena))
  }
}

DNAtoAmino_V2 <- function(cadena) {
  # Crear el data frame de codones
  codones_df <- data.frame(
    codon = c(
      "TTT", "TTC", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG",
      "ATT", "ATC", "ATA", "ATG", "GTT", "GTC", "GTA", "GTG",
      "TCT", "TCC", "TCA", "TCG", "CCT", "CCC", "CCA", "CCG",
      "ACT", "ACC", "ACA", "ACG", "GCT", "GCC", "GCA", "GCG",
      "TAT", "TAC", "TAA", "TAG", "CAT", "CAC", "CAA", "CAG",
      "AAT", "AAC", "AAA", "AAG", "GAT", "GAC", "GAA", "GAG",
      "TGT", "TGC", "TGA", "TGG", "CGT", "CGC", "CGA", "CGG",
      "AGT", "AGC", "AGA", "AGG", "GGT", "GGC", "GGA", "GGG"),
    aminoacid = c(
      "F", "F", "L", "L", "L", "L", "L", "L",
      "I", "I", "I", "M", "V", "V", "V", "V",
      "S", "S", "S", "S", "P", "P", "P", "P",
      "T", "T", "T", "T", "A", "A", "A", "A",
      "Y", "Y", "-", "-", "H", "H", "Q", "Q",
      "N", "N", "K", "K", "D", "D", "E", "E",
      "C", "C", "-", "W", "R", "R", "R", "R",
      "S", "S", "R", "R", "G", "G", "G", "G")
  )
  
  # Función interna para traducir una sola cadena
  traducirCadena <- function(cadena) {
    # Manejar NA en la entrada
    if (is.na(cadena)) {
      return(NA)
    }
    
    if (nchar(cadena) < 1) {
      return("-")  # Para cadenas vacías
    }
    
    aminoacidos <- c()
    
    for (i in seq(1, nchar(cadena), 3)) {
      triplete <- substr(cadena, i, i + 2)
      aa <- codones_df$aminoacid[which(codones_df$codon == triplete)]
      
      if (length(aa) == 0 || nchar(triplete) < 3) {
        aa <- "-"  # Asignar "-" para tripletes no encontrados o incompletos
      }
      
      aminoacidos <- c(aminoacidos, aa)
    }
    
    return(paste(aminoacidos, collapse = ""))
  }
  
  # Aplicar 'traducirCadena' a cada elemento si 'cadena' es un vector
  if(is.vector(cadena)) {
    res <- sapply(cadena, traducirCadena, USE.NAMES = FALSE)
    return(res)
  } else {
    return(traducirCadena(cadena))
  }
}




#--------------------------------------------------#
#         Codon stop                             #
#--------------------------------------------------#
library(dplyr)

# Función para detectar codones de parada
StopCodonHunter <- function(columna) {
  # Codones de parada
  codones_parada <- c("TAA", "TAG", "TGA")
  
  # Aplicar la función a la columna y devolver TRUE si se encuentra algún codón de parada
  sapply(columna, function(secuencia) {
    any(sapply(codones_parada, function(codon) grepl(codon, secuencia)))
  })
}


# Aplicar la función al dataframe
 #df$tiene_codon_parada <- detectar_codon_parada(df$secuencia_nucleotidos)



#--------------------------------------------------#
#           TRaducir de RNA to amino               #
#--------------------------------------------------#

# No aceptar numeros
#Mensaje: No cadena en multiplo de 3

RNAtoAmino <- function(cadena) {
  
  # Crear el data frame de codones
  codones_df <- data.frame(
    codon = c(
      "UUU", "UUC", "UUA", "UUG", "CUU", "CUC", "CUA", "CUG",
      "AUU", "AUC", "AUA", "AUG", "GUU", "GUC", "GUA", "GUG",
      "UCU", "UCC", "UCA", "UCG", "CCU", "CCC", "CCA", "CCG",
      "ACU", "ACC", "ACA", "ACG", "GCU", "GCC", "GCA", "GCG",
      "UAU", "UAC", "UAA", "UAG", "CAU", "CAC", "CAA", "CAG",
      "AAU", "AAC", "AAA", "AAG", "GAU", "GAC", "GAA", "GAG",
      "UGU", "UGC", "UGA", "UGG", "CGU", "CGC", "CGA", "CGG", "AGU", "AGC", "AGA", "AGG",
      "GGU", "GGC", "GGA", "GGG"
    ),
    aminoacid = c(
      "F", "F", "L", "L", "L", "L", "L", "L",
      "I", "I", "I", "M", "V", "V", "V", "V",
      "S", "S", "S", "S", "P", "P", "P", "P",
      "U", "U", "U", "U", "A", "A", "A", "A",
      "Y", "Y", "-", "-", "H", "H", "Q", "Q",
      "N", "N", "K", "K", "D", "D", "E", "E",
      "C", "C", "-", "W", "R", "R", "R", "R", "S", "S", "R", "R",
      "G", "G", "G", "G"
    )
  )
  
  # Función interna para traducir una sola cadena
  traducirCadena <- function(cadena) {
    if(nchar(cadena) < 1) {
      return("-") # o return("") si prefieres una cadena vacía en vez de "-"
    }
    
    aminoacidos <- c()
    
    for (i in seq(1, nchar(cadena), 3)) {
      triplete <- substr(cadena, i, i + 2)
      
      # Buscar el triplete en el data frame de codones
      aa <- codones_df$aminoacid[which(codones_df$codon == triplete)]
      
      # Si el triplete no se encontró o no tiene 3 bases, se asigna "-"
      if (length(aa) == 0 || nchar(triplete) < 3) {
        aa <- "-"
      }
      
      aminoacidos <- c(aminoacidos, aa)
    }
    
    return(paste(aminoacidos, collapse = ""))
  }
  
  # Si 'cadena' es un vector, aplicar 'traducirCadena' a cada elemento
  if(is.vector(cadena)) {
    res <- sapply(cadena, traducirCadena)
    names(res) <- NULL  # quitar nombres del vector
    return(res)
  } else {
    return(traducirCadena(cadena))
  }
}

RNAtoAmino("UUUUUUUGG")


DNAtoAmino <- function(cadena) {
  # Diccionario
  codones_df <- data.frame(
    codon = c(
      "TTT", "TTC", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG",
      "ATT", "ATC", "ATA", "ATG", "GTT", "GTC", "GTA", "GTG",
      "TCT", "TCC", "TCA", "TCG", "CCT", "CCC", "CCA", "CCG",
      "ACT", "ACC", "ACA", "ACG", "GCT", "GCC", "GCA", "GCG",
      "TAT", "TAC", "TAA", "TAG", "CAT", "CAC", "CAA", "CAG",
      "AAT", "AAC", "AAA", "AAG", "GAT", "GAC", "GAA", "GAG",
      "TGT", "TGC", "TGA", "TGG", "CGT", "CGC", "CGA", "CGG",
      "AGT", "AGC", "AGA", "AGG", "GGT", "GGC", "GGA", "GGG"),
    aminoacid = c(
      "F", "F", "L", "L", "L", "L", "L", "L",
      "I", "I", "I", "M", "V", "V", "V", "V",
      "S", "S", "S", "S", "P", "P", "P", "P",
      "T", "T", "T", "T", "A", "A", "A", "A",
      "Y", "Y", "-", "-", "H", "H", "Q", "Q",
      "N", "N", "K", "K", "D", "D", "E", "E",
      "C", "C", "-", "W", "R", "R", "R", "R",
      "S", "S", "R", "R", "G", "G", "G", "G")
  )
  
  # Función interna para traducir una sola cadena
  traducirCadena <- function(cadena) {
    # Manejar NA en la entrada
    if (is.na(cadena)) {
      return(NA)
    }
    
    if (nchar(cadena) < 1) {
      return("-")  # Para cadenas vacías
    }
    
    aminoacidos <- c()
    
    for (i in seq(1, nchar(cadena), 3)) {
      triplete <- substr(cadena, i, i + 2)
      aa <- codones_df$aminoacid[which(codones_df$codon == triplete)]
      
      if (length(aa) == 0 || nchar(triplete) < 3) {
        aa <- "-"  # Asignar "-" para tripletes no encontrados o incompletos
      }
      
      aminoacidos <- c(aminoacidos, aa)
    }
    
    return(paste(aminoacidos, collapse = ""))
  }
  
  # Aplicar 'traducirCadena' a cada elemento si 'cadena' es un vector
  if(is.vector(cadena)) {
    res <- sapply(cadena, traducirCadena, USE.NAMES = FALSE)
    return(res)
  } else {
    return(traducirCadena(cadena))
  }
}


#------------------------------#
#     a la caza de K-MERS      #
#------------------------------#

KmerHunter_V1 <- function(secuencias, k) {
  sapply(secuencias, function(secuencia) {
    n <- nchar(secuencia)
    if (n < k) return("")  
    
    kmers <- character(n - k + 1)
    
    for(i in 1:(n - k + 1)) {
      kmers[i] <- substr(secuencia, i, i + k - 1)
    }
    
    return(paste(kmers, collapse = ";"))
  })
}


KmerHunter <- function(secuencias, k) {
  sapply(secuencias, function(secuencia) {
    # Verificar si la secuencia es NA
    if (is.na(secuencia)) {
      return(NA)  # Devuelve NA directamente si la secuencia es NA
    }
    
    n <- nchar(secuencia)
    if (n < k) return("")  # Retorna una cadena vacía si la secuencia es más corta que k
    
    kmers <- character(n - k + 1)
    
    for(i in 1:(n - k + 1)) {
      kmers[i] <- substr(secuencia, i, i + k - 1)
    }
    
    return(paste(kmers, collapse = ";"))
  })
}

#Añadir errores: Argumento K es missing


#----------------------------------------------------#
#     Extraer sequencias de nucleotidos a fasta      #
#----------------------------------------------------#

library(Biostrings)

# Definición de la función
seqDNAtoFasta <- function(dataframe, columna_nombres, columna_secuencias, ruta_salida) {
  # Extraer secuencias y nombres
  secuencias <- dataframe[[columna_secuencias]]
  nombres <- dataframe[[columna_nombres]]
  
  # Configurar nombres
  names(secuencias) <- nombres
  
  # Convertir a DNAStringSet
  dna <- DNAStringSet(secuencias)
  
  # Escribir al archivo
  writeXStringSet(dna, format = "fasta", file = ruta_salida)
}


#----------------------------------------------------#
#     Extraer sequencias de nucleotidos a fasta      #
#----------------------------------------------------#
# dependencia: Biostrins

seqAAtoFasta <- function(dataframe, columna_nombres, columna_secuencias, ruta_salida) {
  # Extraer secuencias y nombres
  secuencias <- dataframe[[columna_secuencias]]
  nombres <- dataframe[[columna_nombres]]
  
  # Configurar nombres
  names(secuencias) <- nombres
  
  # Convertir a AAStringSet
  aa <- AAStringSet(secuencias)
  
  # Escribir al archivo
  writeXStringSet(aa, format = "fasta", file = ruta_salida)
}

# Uso de la función
# seqAAtoFasta(tu_dataframe, "nombre_columna_nombres", "nombre_columna_secuencias", "ruta_salida.fasta")



#--------------------------------------------------#
#        Extracting VDJ Kappa                      #
#--------------------------------------------------#

DNA_kappaVDJ.Extract <- function(sequence) {
  
  amino_to_codons <- list(
    D = c("GAT", "GAC"),
    I = c("ATT", "ATC", "ATA"),
    V = c("GTT", "GTC", "GTA", "GTG"),
    E = c("GAA", "GAG"),
    L = c("TTA", "TTG", "CTT", "CTC", "CTA", "CTG"),
    K = c("AAA", "AAG")
  )
  
  generate_codon_combinations <- function(amino_sequence) {
    codon_lists <- lapply(strsplit(amino_sequence, NULL)[[1]], function(aa) amino_to_codons[[aa]])
    combinations <- expand.grid(codon_lists)
    return(apply(combinations, 1, paste0, collapse = ""))
  }
  
  kappa_starts <- c(generate_codon_combinations("DIV"), generate_codon_combinations("EIV"))
  kappa_end <- generate_codon_combinations("LEIK")
  
  for(start in kappa_starts) {
    for(end in kappa_end) {
      start_idx <- regexpr(start, sequence)[[1]]
      if(start_idx != -1) {
        end_search_space <- substr(sequence, start_idx + nchar(start), nchar(sequence))
        end_idx <- regexpr(end, end_search_space)[[1]]
        if(end_idx != -1) {
          return(substr(sequence, start_idx, start_idx + nchar(start) + end_idx + nchar(end) - 2))
        }
      }
    }
  }
  
  return(NULL)
}

# Test


Acute_Bcells_dataset_IGL %>% filter(constant == "IGL")
#DNA_kappaVDJExtract(Acute_Bcells_dataset_IGL$sequence)



#-------------------------------------#
#        Extraer lambda
#-------------------------------------#
DNA_lambdaVDJ.Extract <- function(sequence) {
  
  amino_to_codons <- list(
    Q = c("CAA", "CAG"),
    A = c("GCT", "GCC", "GCA", "GCG"),
    V = c("GTT", "GTC", "GTA", "GTG"),
    L = c("TTA", "TTG", "CTT", "CTC", "CTA", "CTG"),
    T = c("ACT", "ACC", "ACA", "ACG")
  )
  
  generate_codon_combinations <- function(amino_sequence) {
    codon_lists <- lapply(strsplit(amino_sequence, NULL)[[1]], function(aa) amino_to_codons[[aa]])
    combinations <- expand.grid(codon_lists)
    return(apply(combinations, 1, paste0, collapse = ""))
  }
  
  lambda_start <- generate_codon_combinations("QAV")
  lambda_end <- generate_codon_combinations("LTVL")
  
  for(start in lambda_start) {
    for(end in lambda_end) {
      start_idx <- regexpr(start, sequence)[[1]]
      if(start_idx != -1) {
        end_search_space <- substr(sequence, start_idx + nchar(start), nchar(sequence))
        end_idx <- regexpr(end, end_search_space)[[1]]
        if(end_idx != -1) {
          return(substr(sequence, start_idx, start_idx + nchar(start) + end_idx + nchar(end) - 2))
        }
      }
    }
  }
  
  return(NULL)
}

# Test
#seq <- "asdfgteerpkbmCAAGCAGTGTACCAGGCCGTTGCGTAGCTTACTGTTCTGTGCTGACTGTCTTATGLSSDFPǴBMSBP"
#print(DNA_lambdaVDJ.Extract(seq))


#------------------------------------------
#          Extraer VDJ Heavy       #
#-----------------------------------#

DNA_HeavyVDJ.extract <- function(df, columna) {
  
  # Definición de los codones
  codones <- list(
    E = c("GAA", "GAG"),
    V = c("GTT", "GTC", "GTA", "GTG"),
    Q = c("CAA", "CAG"),
    L = c("CTT", "CTC", "CTA", "CTG", "TTA", "TTG"),
    T = c("ACT", "ACC", "ACA", "ACG"),
    S = c("TCT", "TCC", "TCA", "TCG", "AGT", "AGC")
  )
  
  # Función para generar todas las combinaciones posibles de una secuencia de aminoácidos
  generar_combinaciones <- function(seq_aminoacidos) {
    combinaciones <- expand.grid(lapply(strsplit(seq_aminoacidos, NULL)[[1]], function(aa) { codones[[aa]] }))
    combinaciones <- apply(combinaciones, 1, paste0, collapse = "")
    return(combinaciones)
  }
  
  # Generar todas las combinaciones posibles para EVQL y TVSS
  combinaciones_EVQL <- generar_combinaciones("EVQL")
  combinaciones_TVSS <- generar_combinaciones("TVSS")
  
  # Inicializar un vector para almacenar las subcadenas encontradas
  subcadenas_encontradas <- vector("character", nrow(df))
  
  # Iterar sobre cada fila en el DataFrame
  for (i in seq_len(nrow(df))) {
    secuencia <- df[i, columna]
    subcadenas_i <- c()  # Inicializar un vector para almacenar las subcadenas encontradas para esta secuencia
    
    # Buscar cada combinación de EVQL
    for (evql in combinaciones_EVQL) {
      pos_inicio <- regexpr(evql, secuencia)
      if (pos_inicio > 0) {  # Si se encontró EVQL
        
        # Buscar cada combinación de TVSS después de EVQL
        for (tvss in combinaciones_TVSS) {
          pos_fin <- regexpr(tvss, substring(secuencia, pos_inicio + nchar(evql)))
          if (pos_fin > 0) {  # Si se encontró TVSS
            subcadena <- substring(secuencia, pos_inicio, pos_inicio + nchar(evql) + pos_fin + nchar(tvss) - 1)
            subcadenas_i <- c(subcadenas_i, subcadena)
          }
        }
      }
    }
    
    # Concatenar las subcadenas encontradas para esta secuencia en una cadena de texto, separadas por ";"
    subcadenas_encontradas[i] <- paste(subcadenas_i, collapse = ";")
  }
  
  return(subcadenas_encontradas)
}



df_test <-
Acute_Bcells_dataset_IGH %>% filter(codigo == 5904) %>% filter(CDR3Length == 54) %>% filter(celltype == "Plasmablast") %>% filter(VGene_clonalyst == "IGHV3-23")
  
DNA_HeavyVDJ.extract(df_test$sequence)
DNA_HeavyVDJ.extract(df_test, "sequence")

#----------------------------------#
#    EXTRAER VDJ A PARTIR DE AA 2 #
#---------------------------------#
AA_HeavyVDJ.extract <- function(df, columna) {
  
  # Definir las secuencias de aminoácidos de interés
  seq_inicio <- "EVQL"
  seq_fin <- "TVSS"
  
  # Inicializar un vector para almacenar las subcadenas encontradas
  subcadenas_encontradas <- vector("character", nrow(df))
  
  # Iterar sobre cada fila en el DataFrame
  for (i in seq_len(nrow(df))) {
    secuencia <- df[i, columna]
    subcadenas_i <- c()  # Inicializar un vector para almacenar las subcadenas encontradas para esta secuencia
    
    # Buscar la secuencia de inicio
    pos_inicio <- regexpr(seq_inicio, secuencia)
    if (pos_inicio > 0) {  # Si se encontró la secuencia de inicio
      
      # Buscar la secuencia de fin después de la secuencia de inicio
      pos_fin <- regexpr(seq_fin, substring(secuencia, pos_inicio + nchar(seq_inicio) - 1))
      if (pos_fin > 0) {  # Si se encontró la secuencia de fin
        # Ajustar pos_fin para reflejar la posición en la cadena original
        pos_fin <- pos_fin + pos_inicio + nchar(seq_inicio) - 2
        subcadena <- substring(secuencia, pos_inicio, pos_fin + nchar(seq_fin) - 1)
        subcadenas_i <- c(subcadenas_i, subcadena)
      }
    }
    
    # Concatenar las subcadenas encontradas para esta secuencia en una cadena de texto, separadas por ";"
    subcadenas_encontradas[i] <- paste(subcadenas_i, collapse = ";")
  }
  
  return(subcadenas_encontradas)
}

          # AA_HeavyVDJ.extract(df, "column")

generateORICP


#-------------------------------------#
#     Extraer fasta facilente        #
#-------------------------------------#
GetFasta <- function(dataframe, header, sequence, output_file = NULL) {
  # Asegurar que los nombres de las columnas existen en el dataframe
  if (!(header %in% names(dataframe)) || !(sequence %in% names(dataframe))) {
    stop("Please ensure the header and sequence column names exist in the dataframe.")
  }
  
  # Crear un vector con las secuencias
  sequences <- Biostrings::DNAStringSet(as.character(dataframe[[sequence]]))
  
  # Establecer los nombres de las secuencias
  names(sequences) <- dataframe[[header]]
  
  # Especificar la ruta del archivo de salida si no se proporciona
  if (is.null(output_file)) {
    output_file <- paste0("output_", Sys.Date(), ".fasta")
  }
  
  # Exportar las secuencias a un archivo FASTA
  writeXStringSet(sequences, filepath = output_file)
  
  # Mensaje de éxito
  cat("FASTA file has been successfully written to", output_file, "\n")
}
