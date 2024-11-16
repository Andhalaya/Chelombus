#!/bin/bash
# Bash script to process data in all files
# Configuración
INPUT_DIR="/mnt/10tb_hdd/uncleaned_data/"
PYTHON_SCRIPT="sort_database.py"
OUTPUT_DIR="/mnt/samsung_2tb/mixed_data"

# Asegúrate de que el script Python existe
if [ ! -f "$PYTHON_SCRIPT" ]; then
  echo "Error: No se encuentra el script $PYTHON_SCRIPT"
  exit 1
fi

# Crear el directorio de salida si no existe
mkdir -p "$OUTPUT_DIR"

# Procesar cada archivo .cxsmiles en el directorio de entrada
for input_file in "$INPUT_DIR"/*.cxsmiles; do
  if [ -f "$input_file" ]; then
    echo "Processing file: $input_file"
    # Ejecutar el script Python para este archivo
    python3 "$PYTHON_SCRIPT" "$input_file"
    if [ $? -ne 0 ]; then
      echo "Error when processing $input_file"
      exit 1
    fi
  else
    echo "No files .cxsmiles files found in $INPUT_DIR"
    exit 1
  fi
done

echo "Processing completed"