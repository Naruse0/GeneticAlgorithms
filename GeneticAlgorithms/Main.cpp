//!	@brief	遺伝的アルゴリズムでナップサック問題を解く

//---- インクルード ----
#include	<stdio.h>
#include	<stdlib.h>
#include	<time.h>

//---- 定数定義 ----
#define		BAG_SIZE	400				// バッグの許容重量
#define		BAG_CAPA	10				// バッグの積載量
#define		ITEM_NUM	64
#define		POP_SIZE	10				// 遺伝子総数
#define		GENE_LENGTH	(ITEM_NUM)		// 遺伝子のビット数

#define		CROS_RATE	0.8f			// 交叉頻度
#define		MUTE_RATE	0.0f			// 突然変異頻度

#define		GENE_MAX	10				// 世代交代数

const float	g_DataVolume[] = { 
  0.072431f, 0.470347f, 0.064719f, 0.285837f, 0.225093f, 0.096393f, 0.515726f, 0.455327f,
  0.165786f, 0.167166f, 0.536347f, 0.147685f, 0.091873f, 0.496223f, 0.253056f, 0.516451f,
  0.067792f, 0.544281f, 0.523754f, 0.277838f, 0.078748f, 0.167894f, 0.302647f, 0.199729f,
  0.398927f, 0.366057f, 0.500194f, 0.192629f,  0.07454f,   0.1858f, 0.288475f, 0.054764f,
  0.178058f, 0.543278f, 0.082055f, 0.147755f, 0.056136f, 0.377898f, 0.149305f, 0.172916f,
  0.547585f, 0.379011f, 0.504197f, 0.078298f, 0.329076f, 0.494272f, 0.528666f, 0.125196f,
  0.515643f, 0.516848f, 0.200623f, 0.470438f, 0.074357f, 0.272055f, 0.421436f, 0.058918f,
  0.530786f, 0.327218f, 0.186787f, 0.488957f, 0.257442f, 0.326924f, 0.143436f, 0.107303f };

const float g_DataWeight[] = { 
  1.042288f, 1.096132f, 0.952216f,  0.79302f, 1.046246f, 0.823767f, 0.873065f, 0.758806f,
  0.737015f, 1.028084f, 0.606415f, 0.881106f, 0.884353f, 1.290026f, 0.774216f, 1.014714f,
  1.229122f, 1.015562f, 0.964802f, 1.288146f, 1.076072f, 0.952694f, 0.704871f, 0.904247f,
  0.846266f, 0.837016f, 1.059204f, 0.699134f, 0.870344f, 0.858246f, 0.743155f, 0.800332f,
  1.160441f, 0.715929f, 1.206835f, 1.056951f, 0.898189f, 1.226303f, 1.025044f, 0.640907f,
  0.899368f, 1.225172f, 0.901907f, 1.068226f, 0.712972f, 1.055511f, 1.114367f, 0.632971f,
  1.235422f, 0.818592f, 1.277626f, 0.976713f, 0.825303f, 0.786353f, 0.904311f, 1.125975f,
	0.9002f, 0.687584f, 1.253057f, 0.609091f,  0.88333f,  0.89934f, 0.928529f, 1.195918f };

const int g_DataValue[] = {
  515, 802, 778, 406, 146, 914, 929, 317,
  556, 788, 438, 897, 127, 425, 909, 650,
  802, 532, 768, 524, 376, 204, 113, 627,
  256, 869, 488, 932, 229, 418, 147, 234,
  210, 290, 161, 655, 820, 981, 835, 903,
  864, 468, 908, 421, 128, 446, 298, 289,
  762, 358, 165, 636, 631, 781, 196, 421,
  446, 402, 431, 205, 511, 441, 680, 688 };


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
	int		i, j, wk_value;
	float	wk_weight, wk_volume;

	// 初期化
	i = 0;
	j = 0;
	g_Generation = 0;
	g_AllMaxFitness = 0;
	g_AllMaxFitnessGeneration = 0;
	g_AllMaxFitnessIndex = 0;

	initialize_pop_binary();

	// 初期状態を表示
	printf("---- 初期状態 ----\n");
	check_result();
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
	printf("最高適応度       = %3d\n", g_AllMaxFitness);
	printf("最高遺伝子世代   = %3d\n", g_AllMaxFitnessGeneration);
	printf("最高遺伝子[%2d]   = ", g_AllMaxFitnessIndex);

	wk_value = 0;
	wk_weight = 0;
	wk_volume = 0;
	for (i = 0; i < GENE_LENGTH; i++)
	{
		printf("%1d", g_AllMaxFitnessGene[i]);
		if (g_AllMaxFitnessGene[i] == 1)
		{
			wk_value += g_DataValue[i];
			wk_weight += g_DataWeight[i];
			wk_volume += g_DataVolume[i];
		}
	}
	printf("\n");
	printf("\n");
	printf("価値合計 = %d\n", wk_value);
	printf("総重量   = %f\n", wk_weight);
	printf("総体積   = %f\n", wk_volume);

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
				g_Gene[index1][j] = gene1[j];
				g_Gene[index2][j] = gene2[j];
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
				//g_Gene[i][pos] != 0 ? g_Gene[i][pos] = 0 : g_Gene[i][pos] = 1;
				g_Gene[i][pos] = !g_Gene[i][pos];
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
			printf("%d", g_AllMaxFitnessGene[i]);
		}
		printf("\n");
		printf("適応度 = %d\n", g_AllMaxFitness);
	}
}

int calc_fitness_gene(int gene_no)
{
	int i, work_value;
	float work_weight, work_volume;

	work_weight = 0;
	work_value = 0;
	work_volume = 0;

	for (i = 0; i < GENE_LENGTH; i++)
	{
		// 遺伝子情報が「１」なら「重さ」「価値」を加算する
		if (g_Gene[gene_no][i] == 1)
		{
			work_weight += g_DataWeight[i];
			work_volume += g_DataVolume[i];
			work_value += g_DataValue[i];
		}
	}

	// もし合計値がバッグの許容量をオーバーしたら適応度０とする
	if (work_weight > BAG_SIZE || 
		work_volume > BAG_CAPA)
	{
		work_value = 0;
	}
	return work_value;
}

void check_result()
{
	int i, j;

	for (i = 0; i < POP_SIZE; i++)
	{
		printf("遺伝子[%2d] : ", i);
		
		for (j = 0; j < GENE_LENGTH; j++)
		{
			printf("%1d", g_Gene[i][j]);
		}
		printf("\n");
	}
}