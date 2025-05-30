# Makefile to compile cavity.c

# コンパイラの指定
CC = gcc

# コンパイルフラグ
# CFLAGS = -Wall -O2
CFLAGS = -O2

# 実行ファイル名
TARGET = cavity

# ソースファイル
SRCS = cavity.c myIO.c myConst.c

# オブジェクトファイル
OBJS = $(SRCS:.c=.o)

# デフォルトのターゲット
all: $(TARGET)

# post process
postProcess: myPostProcess.o myIO.o myConst.o
	$(CC) $(CFLAGS) -o $@ $^

# test
test: test.o myIO.o myConst.o
	$(CC) $(CFLAGS) -o $@ $^

# couette
couette: cavity.o myIO.o myConst.o
	$(CC) $(CFLAGS) -o $@ $^

# poiseuille
poi: cavity.o myIO.o myConst.o
	$(CC) $(CFLAGS) -o $@ $^

# 実行ファイルの生成ルール
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^

# オブジェクトファイルの生成ルール
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# クリーンアップ（オブジェクトファイルと実行ファイルを削除）
clean:
	rm -f *.o $(TARGET)

# cavity.cを再コンパイル
rebuild: clean all
