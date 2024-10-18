# Makefile to compile cavity.c

# コンパイラの指定
CC = gcc

# コンパイルフラグ
# CFLAGS = -Wall -O2
CFLAGS = -O2

# 実行ファイル名
TARGET = cavity

# デフォルトのターゲット
all: $(TARGET)

# cavity.cをコンパイルして実行ファイルを生成
$(TARGET): cavity.o
	$(CC) $(CFLAGS) -o $(TARGET) cavity.o

# cavity.cをコンパイルしてオブジェクトファイルを生成
cavity.o: cavity.c
	$(CC) $(CFLAGS) -c cavity.c

# クリーンアップ（オブジェクトファイルと実行ファイルを削除）
clean:
	rm -f *.o $(TARGET)

# cavity.cを再コンパイル
rebuild: clean all
