#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="$SCRIPT_DIR/build"
USER_ID=$(id -u)
GROUP_ID=$(id -g)

echo "👤 Running as user: $USER_ID:$GROUP_ID"

docker build \
    --build-arg USER_ID=$USER_ID \
    --build-arg GROUP_ID=$GROUP_ID \
    -f "$SCRIPT_DIR/Dockerfile.build" \
    -t build-package "$SCRIPT_DIR/"

docker run \
    --user "$USER_ID:$GROUP_ID" \
    -v "$BUILD_DIR:/output" \
    build-package

echo "✅ Build complete! Files owned by $(whoami)"

# Копирование собранных файлов в целевые директории
echo "📁 Copying build files to target directories..."

# Определяем целевые пути относительно корневой директории TGC
TGC_ROOT="$SCRIPT_DIR/../.."
FILE_WORKER_WHL="$TGC_ROOT/file-worker/whl"
CALC_WORKER_WHL="$TGC_ROOT/calc-worker/whl"
WEBSITE_WHL="$TGC_ROOT/website/whl"

# Создаем целевые директории, если они не существуют
mkdir -p "$FILE_WORKER_WHL"
mkdir -p "$CALC_WORKER_WHL"
mkdir -p "$WEBSITE_WHL"

# Копируем содержимое build директории
echo "📋 Copying to file-worker/whl..."
cp -r "$BUILD_DIR"/* "$FILE_WORKER_WHL/" 2>/dev/null || true

echo "📋 Copying to calc-worker/whl..."
cp -r "$BUILD_DIR"/* "$CALC_WORKER_WHL/" 2>/dev/null || true

echo "📋 Copying to website/whl..."
cp -r "$BUILD_DIR"/* "$WEBSITE_WHL/" 2>/dev/null || true

echo "🎉 All files copied successfully!"