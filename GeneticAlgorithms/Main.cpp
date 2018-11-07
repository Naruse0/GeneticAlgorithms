//!	@brief	遺伝的アルゴリズムでナップサック問題を解く

//---- インクルード ----
#include	<stdio.h>
#include	<stdlib.h>
#include	<time.h>

//---- 定数定義 ----
#define		BAG_SIZE	400
#define		ITEM_NUM	20
#define		POP_SIZE	10				// 遺伝子総数
#define		GENE_LENGTH	(ITEM_NUM)		// 遺伝子のビット数

#define		CROS_RATE	0.8f			// 交叉頻度
#define		MUTE_RATE	0.1f			// 突然変異頻度

#define		GENE_MAX	10				// 世代交代数


const int g_DataWeight[] = { 80,70,60,34,89,10,20,30,90,38,75,24,55,15,28,40,50,64,15,79 };
const int g_DataValue[] = { 10, 5, 20,8,2,5,7,14,22,9,11,3,18,13,6,21,24,1,4,16 };


//---- グローバル変数 ----
int g_Generation;
int	g_Gene[POP_SIZE][GENE_LENGTH];
int	g_Fitness[POP_SIZE];				// 適応度集団
int g_AllMaxFitness;					// 最高適応度
int g_AllMaxFitnessGeneration;			// 最高適応度出現世代
int g_AllMaxFitnessIndex;				// 最高適応度遺伝子番号
int g_AllMaxFitnessGene[GENE_LENGTH];	// 最高適応度遺伝子（これが解となる）

//---- プロトタイプ宣言 ----
void initialize_pop_binary();			// 初期化
void M_selection();						// 選択
void M_crossover();						// 交叉
void M_mutation();						// 突然変異
void calc_fitness_pop();				// 遺伝子の適応度を計算
void check_result();					// 現在の遺伝子情報を表示
int calc_fitness_gene(int gene_no);

//---- エントリポイント ---- 
int main()
{
	int	i, j, wk_weight;

	// 初期化
	i = 0;
	j = 0;
	g_Generation = 0;
	g_AllMaxFitness = 0;
	g_AllMaxFitnessGeneration = 0;
	g_AllMaxFitnessIndex = 0;

	initialize_pop_binary();

	// 初期状態を表示

	// 遺伝子集団の適応度を計算
	calc_fitness_pop();

	// 世代交代をしながら遺伝子を進化させる
	for (g_Generation = 1; g_Generation < GENE_MAX; g_Generation++)
	{
		printf("---- 第 %3d 世代 ----\n", g_Generation);
		M_selection();
		M_crossover();
		M_mutation();
		calc_fitness_pop();
		check_result();
	}

	// 結果を出力
	printf("==== 計算終了 ====\n");
	printf("最高適応度		= %3d\n", g_AllMaxFitness);
	printf("最高遺伝子世代	= %3d\n", g_AllMaxFitnessGeneration);
	printf("最高遺伝子[%2d]	= ", g_AllMaxFitnessIndex);

	wk_weight = 0;
	for (i = 0; i < GENE_LENGTH; i++)
	{
		printf(" %1d", g_AllMaxFitnessGene[i]);
		if (g_AllMaxFitnessGene[i] == 1)
		{
			wk_weight += g_DataWeight[i];
		}
	}
	printf("\n");
	printf("総重量 = %d\n", wk_weight);


	// 終了処理
	printf("Push Enter-key >> ");
	rewind(stdin);
	getchar();
	return 0;
}

//---- 関数 ----
void initialize_pop_binary()
{
	int i, j;

	srand(static_cast<unsigned int>(time(nullptr)));
	
	for (i = 0; i < POP_SIZE; i++)
	{
		for (j = 0; j < GENE_LENGTH; j++)
		{
			// 0,1で初期化
			int geneValue = 0;
			rand() - static_cast<int>(RAND_MAX * 0.5f) >= 0 ? geneValue = 1 : geneValue = 0;
			g_Gene[i][j] = geneValue;
		}
	}
}

void M_selection()
{
	int i, j;

	// 2つずつ適応度を比較し、適応度の高いほうの遺伝子をもう片方にコピーする
	for (i = 0; i < (POP_SIZE - 1); i += 2)
	{
		if (g_Fitness[i] < g_Fitness[i + 1])
		{
			for (j = 0; j < GENE_LENGTH; j++)
			{
				g_Gene[i][j] = g_Gene[i + 1][j];
			}
		}
	}
}

void M_crossover()
{
	int				gene1[GENE_LENGTH], gene2[GENE_LENGTH];
	unsigned char	work;
	int				i, j;
	int				c_pos;
	double			r;

	int				indices[POP_SIZE];	// 処理待ちのインデックス配列(-1で処理済み)

	// インデックス配列を初期化
	for (i = 0; i < POP_SIZE; i++)
	{
		indices[i] = i;
	}


	for (i = 0; i < (POP_SIZE - 1); i += 2)
	{
		// 交叉させるかどうかを乱数から判定する
		r = (static_cast<double>(rand() % 10001) / 10000.0);
		if (r <= CROS_RATE)
		{
			// 処理対象のインデックスを決める
			int index1;
			int index2;
			do
			{
				index1 = indices[rand() % POP_SIZE];
			} while (index1 == -1);
			do
			{
				index2 = indices[rand() % POP_SIZE];
			} while (index2 == -1 || index2 == index1);

			// 処理済みとしてカウント
			indices[index1] = -1;
			indices[index2] = -1;

			// 一度に交叉させる2つの遺伝子をワークにコピーする
			for (j = 0; j < GENE_LENGTH; j++)
			{
				gene1[j] = g_Gene[index1][j];
				gene2[j] = g_Gene[index2][j];
			}

			// 乱数を用いて交叉位置を決定し、その値をc_posへ代入する
			c_pos = rand() % GENE_LENGTH;

			// 値以降の遺伝子情報を交叉させる
			for (j = c_pos; j < GENE_LENGTH; ++j)
			{
				work = gene1[j];
				gene1[j] = gene2[j];
				gene2[j] = work;
			}


			// 交叉した新しい遺伝子を遺伝子情報に上書きコピー
			for (j = 0; j < GENE_LENGTH; j++)
			{
				g_Gene[index1][j] = gene1[i];
				g_Gene[index2][j] = gene2[i];
			}
		}
	}
}

void M_mutation()
{
	int		i, j;
	double	r;
	int		pos;

	for (i = 0; i < POP_SIZE; i++)
	{
		for (j = 0; j < GENE_LENGTH; j++)
		{
			// 突然変異が起きるかを乱数で決定
			r = (static_cast<double>(rand() % 10001) / 10000);

			if (r < MUTE_RATE)
			{
				pos = rand() % GENE_LENGTH;

				// 数値を反転（突然変異）
				g_Gene[i][j] == 0 ? g_Gene[i][j] = 1 : g_Gene[i][j] = 0;
			}
		}
	}
}

void calc_fitness_pop()
{
	int	i, work_i, index;
	work_i = 0;
	index = 0;

	for (i = 0; i < POP_SIZE; i++)
	{
		// 指定遺伝子の適応度を計算
		g_Fitness[i] = calc_fitness_gene(i);

		// 適応度の最大値を求める
		if (g_Fitness[i] > work_i)
		{
			work_i = g_Fitness[i];
			index = i;
		}
	}

	// 今までの最高適応度よりも上回ったら更新
	if (work_i > g_AllMaxFitness)
	{
		 // 最高遺伝子情報を保存
		g_AllMaxFitness = work_i;
		g_AllMaxFitnessGeneration = g_Generation;
		g_AllMaxFitnessIndex = index;

		for (i = 0; i < GENE_LENGTH; i++)
		{
			g_AllMaxFitnessGene[i] = g_Gene[g_AllMaxFitnessIndex][i];
		}
		
		// 情報を出力
		printf("良い遺伝子を発見。\n");
		printf("遺伝子[%2d] = ", g_AllMaxFitnessIndex);
		for (i = 0; i < GENE_LENGTH; ++i)
		{
			printf("%2d", g_AllMaxFitnessGene[i]);
		}
		printf("\n");
		printf("適応度 = %d\n", g_AllMaxFitness);
	}
}

int calc_fitness_gene(int gene_no)
{
	int i, work_w, work_v;
	work_w = 0;
	work_v = 0;

	for (i = 0; i < GENE_LENGTH; i++)
	{
		// 遺伝子情報が「１」なら「重さ」「価値」を加算する
		if (g_Gene[gene_no][i] == 1)
		{
			work_w += g_DataWeight[i];
			work_v += g_DataValue[i];
		}
	}

	// もし合計値がバッグの許容量をオーバーしたら適応度０とする
	if (work_w > BAG_SIZE)
	{
		work_v = 0;
	}
	return work_v;
}

void check_result()
{
	int i, j;

	for (i = 0; i < POP_SIZE; i++)
	{
		printf("遺伝子[%2d] : ", i);
		
		for (j = 0; j < GENE_LENGTH; j++)
		{
			printf(" %1d", g_Gene[i][j]);
		}
		printf("\n");
	}
}