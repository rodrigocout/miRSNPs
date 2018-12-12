# ::Objetivos da minha pesquisa::
# Existem dezenas de SNPs - tagSNPs, que marcam locis distintos em doenças 
# complexas. No entanto, a relaçao causal entre o aumento de risco genético e
# início/progressão da doença; e os mecanismos que os polimorfismos exercem 
# suas funções, é desconhecido.

## Pontos Importantes: 
# -> como ocorre a distinção entre os SNPs que sao importantes funcionais 
# e os nao funcionais. [...]

# -> tentar priorizar características genomicas dos SNPs funcionais.

# -> Biofeatures: sao regioes de importancia funcional no genoma, ex.: CDS, 
# mRNA, rRNA, tRNA, regioes de transcrição (TSS), etc.

#----------------------------------------------------------------------------

##Retorna um caminho absoluto que representa o diretório
#de trabalho atual do processo de R
getwd()

## É usado para definir o diretório de trabalho para dir
setwd("/dados/ResearchProjects/mikely/dados_cancer_ovario/")

##ovarian <- file.path("/home/mikely")

# lista objetos do meu R
ls()
[1] "TSS.human.GRCh37" "ovarian"          "ovarian.anno"     "ovarian.bio"     
[5] "ovarian.snp"      "txdb" 

# carregando o pacote FunciSNP
library("FunciSNP")

head(ovarian)
head(ovarian.anno)
head(ovarian.snp)
summary(ovarian)
# armazena o ARQUIVO TABELA2-CANCER-OVARIAN.txt em ovarian.snp 
# armazenando o caminho do arquivo "TABELA2-CANCER-OVARIAN.txt", que contém as informações de cromossomos relacionadas às SNPs-tagSNP e 
# a etnia da população que mais apresentam essas características
ovarian.snp <- read.delim("/dados/ResearchProjects/mikely/dados_cancer_ovario/TABELA2-CANCER-OVARIAN.txt", quote = ""
                          header = F) # Digita o sep para separar os espaços
                                              # da coluna e colocar o nome das colunas
head(ovarian.snp)
## print::
"  V1         V2  V3
1 2:177043226  rs6755777 EUR
2 3:156435952 rs62274042 EUR
3   5:1279790 rs10069690 EUR
4  6:28434693  rs2191035 EUR
5  8:82668568 rs76837345 EUR
6 8:129543949 rs10088218 EUR"


ovarian.snp  ## armazena arquivo com colunas mais espaçadas e organizadas
## print::
"   V1          V2  V3
1  2:177043226   rs6755777 EUR
2  3:156435952  rs62274042 EUR
3    5:1279790  rs10069690 EUR
4   6:28434693   rs2191035 EUR
5   8:82668568  rs76837345 EUR
6  8:129543949  rs10088218 EUR
7   9:16913036   rs4631563 EUR
8  9:136141620   rs2519093 EUR
9  10:21827796   rs1802669 EUR
10 17:36096515    rs757210 EUR
11 17:43569659 rs146746174 EUR
12 17:46500673   rs7207826 EUR
13 19:17390291   rs4808075 EUR
14  1:38072798  rs61776206 EUR
15 4:119947188   rs7671665 EUR
16 4:165908721   rs4691139 EUR
17 17:43516402  rs17631303 EUR
18 8:129601486   rs2165805 EUR
19 17:36099840  rs11651755 EUR"


getwd()
# print:
"home/mikely"

# armazena o CAMINHO da pasta onde guardara o arquivo a ser estudado
# na pasta "peaks" contém os biofeatures das células de câncer de ovário, arquivos com formato
# em extensão ".bed", que são as linhagens celulares do câncer de ovário
# imprime os arquivos de extensão ".bed"
# esses arquivos representam cada biofeature, que sao armazenados 
# numa pasta ("peaks") contendo a extensao ".bed". Esses arquivos que possuem
# algumas caracteristicas importantes: localização de pontos de recursos 
# biologicos, e sao separados em formato ".bed".
ovarian.bio <- ("/dados/ResearchProjects/mikely/nice_name/peaks/") ## se nao der certo coloca-se isso: "/dados/ResearchProjects/mikely/nice_name/peaks/"
ovarian.bio
## print::
"/dados/ResearchProjects/mikely/nice_name/peaks/"


# Este comando lista todos os arquivos que possuem sobreposiçoes com cada
# biofeature de interesse do usuario, i. e., com a extensão a ".bed"
# este comando é para listar os arquivos que estão no diretório ovarian.bio
as.matrix(list.files(ovarian.bio, pattern=".bed$"))
## print::
"[,1]                                      
[1,] "iEEC16_FAIRE_NA_rep1_1.filtered.bed"     
[2,] "iEEC16_FAIRE_NA_rep2_1.filtered.bed"     
[3,] "iEEC16_H3K27ac_NA_rep1_1.filtered.bed"   
[4,] "iEEC16_H3K27ac_NA_rep2_1.filtered.bed"   
[5,] "iEEC16_H3K4me1_NA_rep1_1.filtered.bed"   
[6,] "iEEC16_H3K4me1_NA_rep2_1.filtered.bed"   
[7,] "iFTSEC246_FAIRE_NA_rep1_1.filtered.bed"  
[8,] "iFTSEC246_FAIRE_NA_rep2_1.filtered.bed"  
[9,] "iFTSEC246_H3K27ac_NA_rep1_1.filtered.bed"
[10,] "iFTSEC246_H3K27ac_NA_rep2_1.filtered.bed"
[11,] "iFTSEC246_H3K4me1_NA_rep1_1.filtered.bed"
[12,] "iFTSEC246_H3K4me1_NA_rep2_1.filtered.bed"
[13,] "iFTSEC33_FAIRE_NA_rep1_1.filtered.bed"   
[14,] "iFTSEC33_FAIRE_NA_rep2_1.filtered.bed"   
[15,] "iFTSEC33_H3K27ac_NA_rep1_1.filtered.bed" 
[16,] "iFTSEC33_H3K27ac_NA_rep2_1.filtered.bed" 
[17,] "iFTSEC33_H3K4me1_NA_rep1_1.filtered.bed" 
[18,] "iFTSEC33_H3K4me1_NA_rep2_1.filtered.bed" 
[19,] "iOSE11_FAIRE_NA_rep1_1.filtered.bed"     
[20,] "iOSE11_FAIRE_NA_rep2_1.filtered.bed"     
[21,] "iOSE11_H3K27ac_NA_rep1_1.filtered.bed"   
[22,] "iOSE11_H3K27ac_NA_rep2_1.filtered.bed"   
[23,] "iOSE11_H3K4me1_NA_rep1_1.filtered.bed"   
[24,] "iOSE11_H3K4me1_NA_rep2_1.filtered.bed"   
[25,] "iOSE4_FAIRE_NA_rep1_1.filtered.bed"      
[26,] "iOSE4_FAIRE_NA_rep2_1.filtered.bed"      
[27,] "iOSE4_H3K27ac_NA_rep1_1.filtered.bed"    
[28,] "iOSE4_H3K27ac_NA_rep2_1.filtered.bed"    
[29,] "iOSE4_H3K4me1_NA_rep1_1.filtered.bed"    
[30,] "iOSE4_H3K4me1_NA_rep2_1.filtered.bed"


# ANTES ESTAVA ARMAZENANDO O PRÓPRIO ARQUIVO NESTA VARIAVEL, MAS 
# AGORA ESTÁ ARMAZENANDO O CAMINHO da pasta onde estarão os SNPs
# armazena o CAMINHO DA PASTA onde esta o arquivo, PORQUE DEPOIS É QUE PRECISARÁ DESTA VARIÁVEL PARA ATUAR NO CÓDIGO, NA PARTE DE
# DE INTEGRAÇÃO DOS DADOS os SNPs com os beofeatures
# TABELA2-CANCER-OVARIAN.txt
ovarian.snp <- ("/dados/ResearchProjects/mikely/TABELA2-CANCER-OVARIAN.txt") #==> MUITO IMPORTANTE!!! NAO COLOCAR "/" NO FINAL DE ARQUIVO PQ SENAO NAO IRÁ RECONHECER

# onde está salvo
## print::

ovarian.snp
[1] "/dados/ResearchProjects/mikely/TABELA2-CANCER-OVARIAN.txt"


library(FunciSNP)

# Ou seja, em ovarian.snp contêm arquivos com as 
# seguintes informações dos snp's nos cromossomos:
# associado com um determinado fenótipo
# (por exemplo, doença, característica),
# [1]=> posição do cromossomo, seguida da posição do snp neste cromossomo;
# [2]=> referencia de um tagSNP para este snp "identidade do snp";
# [3]=> etnia da população que mais apresentam as características 
# relacionadas a estes snp's;

# ARQUIVOS :: (".BED" - EXTENSÃO) 
# E em ovarian.bio contêm arquivos com um conjunto de picos - 'peaks', 
# chamado de biofeatures das regiões celulares de câncer de ovário, que
# são regiões enriquecidas epigeneticamente (regiões de picos), 
# sendo que os picos de chip-seq identificados estão muito relacionados
#  com a função biológica, (dados usados para depois fazer anotações de 
# genes - distância mais próxima da TSS, o mais próximo do lincRNA, 
# caracterização genômica, etc...- e anotações de biofeatures - picos de 
# polimorfismos). 



# Adicionando em ovarian informações de cromossomos,
# locais de tagSNP, posições, porcentagem e biofeatures.
# Aqui, ocorre a integração das regioes de snp e de biofeatures, armazenando
# na variavel ovarian. 

# Aqui, cada biofeature é usado para se correlacionar com SNP's, onde os 
# arquivos dos beofeatures devem estar em formato ".bed" e armazenado em uma pasta (peaks);
# alem de outra pasta com arquivos contendo informações relacionadas
# com SNP's da linhagem celular da doença de interesse, neste caso do 
# câncer de ovário.
# RETORNA :: SNPs correlacionados (a partir dos genomas 1000 db) que estão em 
# desequilíbrio de ligação (LD) para uma doença conhecida associada a um tag-SNP
# e sobreposições das características biológicas de cromatina - biofeatures. 
# Estes SNPs correlacionados identificados são caracterizados como SNPs 
# putativos (== suposto) funcionais para uma determinada característica. 

# ::RETORNA A INTEGRAÇÃO DOS TAGSNPS COM OS BEOFEATURES:: CORRELAÇÃO DE SNPS NAS CARACTERISTICAS BIOLOGICAS
# Lista TagSNP com 21 SNPs Tag com 12.532 SNPs nas proximidades, potencialmente correlacionados, sobrepondo, pelo menos, um biofeature
# retorna SNPs correlatas (a partir dos genomas db 1000) que estão em desequilíbrio de ligação (LD) 
# para uma doença conhecida associada tag-SNP e sobreposições características biológicas de cromatina. Estes 
# SNPs correlacionados identificados são caracterizados como SNPs funcionais putativos para um traço particular.		

ovarian <- getFSNPs(snp.regions.file=ovarian.snp, bio.features.loc=ovarian.bio)  # argumentos: caminhos dos arquivos => snps e beofeatures, respectivamente

# print:"
| | _  |  _  _ __  _    _|_ _ 
|^|(/_ | (_ (_)|||(/_    |_(_)

__             __    _ 
|_    __  _  o (_ |\||_)
|  |_|| |(_  | __)| ||  

Version: 1.8.0
System: Linux
::args used::
verbose:                          FALSE
cores in use:                     4
snp.regions.file:                 /dados/ResearchProjects/mikely/TABELA2-CANCER-OVARIAN.txt
p-value adjustment by:            BH
Bio Features:                     35: iEEC16_FAIRE_NA_rep1_1.filtered, iEEC16_FAIRE_NA_rep2_1.filtered, iEEC16_H3K27ac_NA_rep1_1.filtered, iEEC16_H3K27ac_NA_rep2_1.filtered, iEEC16_H3K4me1_NA_rep1_1.filtered, iEEC16_H3K4me1_NA_rep2_1.filtered, iFTSEC246_FAIRE_NA_rep1_1.filtered, iFTSEC246_FAIRE_NA_rep2_1.filtered, iFTSEC246_H3K27ac_NA_rep1_1.filtered, iFTSEC246_H3K27ac_NA_rep2_1.filtered, iFTSEC246_H3K4me1_NA_rep1_1.filtered, iFTSEC246_H3K4me1_NA_rep2_1.filtered, iFTSEC33_FAIRE_NA_rep1_1.filtered, iFTSEC33_FAIRE_NA_rep2_1.filtered, iFTSEC33_H3K27ac_NA_rep1_1.filtered, iFTSEC33_H3K27ac_NA_rep2_1.filtered, iFTSEC33_H3K4me1_NA_rep1_1.filtered, iFTSEC33_H3K4me1_NA_rep2_1.filtered, iOSE11_FAIRE_NA_rep1_1.filtered, iOSE11_FAIRE_NA_rep2_1.filtered, iOSE11_H3K27ac_NA_rep1_1.filtered, iOSE11_H3K27ac_NA_rep2_1.filtered, iOSE11_H3K4me1_NA_rep1_1.filtered, iOSE11_H3K4me1_NA_rep2_1.filtered, iOSE4_FAIRE_NA_rep1_1.filtered, iOSE4_FAIRE_NA_rep2_1.filtered, iOSE4_H3K27ac_NA_rep1_1.filtered, iOSE4_H3K27ac_NA_rep2_1.filtered, iOSE4_H3K4me1_NA_rep1_1.filtered, iOSE4_H3K4me1_NA_rep2_1.filtered, CTCF_only.known, EncodeDnaseI_only.known, EncodeDnaseI_withCTCF.known, EncodeFaire.known, knownGene.Promoters.known
Number of TagSNPs Interrogated:   19 representing 19 unique tagSNPs
Error in grid.Call.graphics(L_upviewport, as.integer(n)) : 
cannot pop the top-level viewport ('grid' and 'graphics' output mixed?)
Error in grid.Call.graphics(L_upviewport, as.integer(n)) : 
cannot pop the top-level viewport ('grid' and 'graphics' output mixed?)
Putative Functional SNPs identified!!
Annotation will begin
~~
Adding lincRNA ... done
Adding gene annotations ... done

Adding genomic annotations ... done

Now do the Funci Dance!"

head(ovarian)
"
TagSNP List with  19  Tag SNPs and 
12532 nearby,  potentially correlated SNPs, that overlap at least one biofeature 
$`R squared: 0.1`
Total R.sq>=0.1 Percent
tagSNPs        19        19  100.00
1K SNPs     12532       765    6.10
Biofeatures    35        33   94.29

$`R squared: 0.5`
Total R.sq>=0.5 Percent
tagSNPs        19        19  100.00
1K SNPs     12532       279    2.23
Biofeatures    35        32   91.43

$`R squared: 0.9`
Total R.sq>=0.9 Percent
tagSNPs        19        13   68.42
1K SNPs     12532       170    1.36
Biofeatures    35        32   91.43
"


# Resumo da TSlist ovarian
# sobreposição de pelo menos x biofeatures, por Tag SNP em um determinado
# (corte) R score $ `quadrado R quadrado: 0,1 em 19 SNPs com Tag um total.
# Mostra na coluna a quantidade de tagSNPs, e na linha superior a quantidade
# de biofeatures sobrepostos por estes tagSNPs.
# Mostra a quantidade de tagSNPs em cada biofeature dos cromossomos listados
# acima, os chamados SNPs recentemente identificados - YAFSNPs

summary(ovarian)   # Resumo da TSlist ovarian
"
TagSNP List with  19  Tag SNPs and 
 12532 nearby,  potentially correlated SNPs, that overlap at least one biofeature 
Number of potentially correlated SNPs 
overlapping at least x biofeatures, per Tag SNP at a specified R squared
$`R squared: 0.1 in 19 Tag SNPs with a total of `
bio.1 bio.2 bio.3 bio.4 bio.5 bio.6 bio.7 bio.8 bio.9 bio.10 bio.11 bio.12
rs10069690          3     2     1     0     0     0     0     0     0      0      0      0
rs10088218         93    45    36    31    28    28    28    26    26     25     23     17
rs11651755         16    11     8     6     5     5     4     2     0      0      0      0
rs146746174        92    65    56    51    43    41    41    39    35     32     31     30
rs17631303         66    39    32    27    19    15    15    13     8      6      5      4
rs1802669          19    15    11    11     9     9     9     9     7      6      4      3
rs2165805          29    11     8     5     3     1     1     0     0      0      0      0
rs2191035          30    23    20    19    16    16    15    13    11     10      7      7
rs2519093          21     7     7     7     6     3     3     3     3      3      1      1
rs4631563          44    35    25    24    19    19    18    18    18     14     10     10
rs4691139          50    28    19    14    14    14    13    13    12     11      9      6
rs4808075          81    67    61    50    45    41    35    33    30     26     22     21
rs61776206         22    13    10     7     7     7     6     6     6      6      4      3
rs62274042         81    73    71    62    54    50    45    38    37     35     31     27
rs6755777          62    50    45    33    29    24    19    17    13     12      9      8
rs7207826          67    49    37    31    24    22    19    18    15     12     10     10
rs757210           13    10     8     7     6     5     4     2     0      0      0      0
rs7671665          33    26    19    15    12    10     9     5     4      4      3      3
rs76837345         12     9     7     6     6     4     3     2     2      1      1      0
TOTAL # 1kgSNPs   834   578   481   406   345   314   287   257   227    203    170    150
bio.13 bio.14 bio.15 bio.16 bio.17 bio.18 bio.19 bio.20 bio.21 bio.22 bio.23
rs10069690           0      0      0      0      0      0      0      0      0      0      0
rs10088218          13     11     11     11     11      8      2      2      1      1      0
rs11651755           0      0      0      0      0      0      0      0      0      0      0
rs146746174         26     19     17     15     13     10     10      7      0      0      0
rs17631303           3      1      0      0      0      0      0      0      0      0      0
rs1802669            3      3      2      1      1      1      1      1      1      1      1
rs2165805            0      0      0      0      0      0      0      0      0      0      0
rs2191035            6      5      3      1      1      0      0      0      0      0      0
rs2519093            1      1      1      0      0      0      0      0      0      0      0
rs4631563            8      7      4      1      1      0      0      0      0      0      0
rs4691139            5      5      5      4      4      4      2      2      2      2      2
rs4808075           19     19     16     14     10      8      5      0      0      0      0
rs61776206           2      2      1      1      1      1      1      1      0      0      0
rs62274042          24     21     18     13     12      8      5      2      2      1      1
rs6755777            7      2      1      0      0      0      0      0      0      0      0
rs7207826            7      6      4      4      3      3      3      2      2      0      0
rs757210             0      0      0      0      0      0      0      0      0      0      0
rs7671665            2      1      1      1      1      1      1      1      1      1      0
rs76837345           0      0      0      0      0      0      0      0      0      0      0
TOTAL # 1kgSNPs    126    103     84     66     58     44     30     18      9      6      4
bio.24 bio.25
rs10069690           0      0
rs10088218           0      0
rs11651755           0      0
rs146746174          0      0
rs17631303           0      0
rs1802669            0      0
rs2165805            0      0
rs2191035            0      0
rs2519093            0      0
rs4631563            0      0
rs4691139            0      0
rs4808075            0      0
rs61776206           0      0
rs62274042           1      1
rs6755777            0      0
rs7207826            0      0
rs757210             0      0
rs7671665            0      0
rs76837345           0      0
TOTAL # 1kgSNPs      1      1

$`R squared: 0.5 in 19 Tag SNPs with a total of `
bio.1 bio.2 bio.3 bio.4 bio.5 bio.6 bio.7 bio.8 bio.9 bio.10 bio.11 bio.12
rs10069690          2     2     1     0     0     0     0     0     0      0      0      0
rs10088218         54    32    25    22    21    21    21    20    20     20     19     14
rs11651755         11     8     7     6     5     5     4     2     0      0      0      0
rs146746174        60    36    30    25    17    15    15    13    10      8      8      7
rs17631303         55    31    25    21    14    12    12    10     7      5      5      4
rs1802669          16    12    11    11     9     9     9     9     7      6      4      3
rs2165805           2     1     0     0     0     0     0     0     0      0      0      0
rs2191035           2     1     1     1     0     0     0     0     0      0      0      0
rs2519093           2     0     0     0     0     0     0     0     0      0      0      0
rs4631563           5     3     1     1     1     1     1     1     1      0      0      0
rs4691139           6     6     5     4     4     4     4     4     3      3      2      2
rs4808075           7     6     6     5     4     4     4     4     3      3      2      2
rs61776206          7     4     2     1     1     1     1     1     1      1      1      1
rs62274042         71    67    65    57    49    46    42    36    35     33     29     25
rs6755777          12    10    10     7     4     3     0     0     0      0      0      0
rs7207826          14     9     5     5     3     3     3     3     3      3      3      3
rs757210            6     4     3     2     2     2     2     2     0      0      0      0
rs7671665           3     3     2     2     1     0     0     0     0      0      0      0
rs76837345          2     1     1     0     0     0     0     0     0      0      0      0
TOTAL # 1kgSNPs   337   236   200   170   135   126   118   105    90     82     73     61
bio.13 bio.14 bio.15 bio.16 bio.17 bio.18 bio.19 bio.20 bio.21 bio.22 bio.23
rs10069690           0      0      0      0      0      0      0      0      0      0      0
rs10088218          11     10     10     10     10      7      2      2      1      1      0
rs11651755           0      0      0      0      0      0      0      0      0      0      0
rs146746174          5      3      2      2      1      0      0      0      0      0      0
rs17631303           3      1      0      0      0      0      0      0      0      0      0
rs1802669            3      3      2      1      1      1      1      1      1      1      1
rs2165805            0      0      0      0      0      0      0      0      0      0      0
rs2191035            0      0      0      0      0      0      0      0      0      0      0
rs2519093            0      0      0      0      0      0      0      0      0      0      0
rs4631563            0      0      0      0      0      0      0      0      0      0      0
rs4691139            2      2      2      2      2      2      2      2      2      2      2
rs4808075            2      2      2      2      2      1      1      0      0      0      0
rs61776206           1      1      0      0      0      0      0      0      0      0      0
rs62274042          22     19     16     11     10      6      4      1      1      1      1
rs6755777            0      0      0      0      0      0      0      0      0      0      0
rs7207826            3      2      0      0      0      0      0      0      0      0      0
rs757210             0      0      0      0      0      0      0      0      0      0      0
rs7671665            0      0      0      0      0      0      0      0      0      0      0
rs76837345           0      0      0      0      0      0      0      0      0      0      0
TOTAL # 1kgSNPs     52     43     34     28     26     17     10      6      5      5      4
bio.24 bio.25
rs10069690           0      0
rs10088218           0      0
rs11651755           0      0
rs146746174          0      0
rs17631303           0      0
rs1802669            0      0
rs2165805            0      0
rs2191035            0      0
rs2519093            0      0
rs4631563            0      0
rs4691139            0      0
rs4808075            0      0
rs61776206           0      0
rs62274042           1      1
rs6755777            0      0
rs7207826            0      0
rs757210             0      0
rs7671665            0      0
rs76837345           0      0
TOTAL # 1kgSNPs      1      1

$`R squared: 0.9 in 13 Tag SNPs with a total of `
bio.1 bio.2 bio.3 bio.4 bio.5 bio.6 bio.7 bio.8 bio.9 bio.10 bio.11 bio.12
rs10088218         31    18    17    16    16    16    16    16    16     16     15     11
rs11651755          4     4     3     2     2     2     2     2     0      0      0      0
rs146746174        40    22    18    15    11    10    10     8     5      3      3      3
rs17631303         46    24    19    16    12    11    11     9     6      4      4      4
rs1802669           6     5     4     4     4     4     4     4     3      2      1      1
rs4631563           1     0     0     0     0     0     0     0     0      0      0      0
rs4691139           1     1     1     1     1     1     1     1     1      1      1      1
rs4808075           5     4     4     3     2     2     2     2     1      1      1      1
rs61776206          1     0     0     0     0     0     0     0     0      0      0      0
rs62274042         57    54    52    45    39    37    35    29    29     28     25     22
rs6755777           9     7     7     5     3     2     0     0     0      0      0      0
rs7207826           7     5     3     3     1     1     1     1     1      1      1      1
rs76837345          2     1     1     0     0     0     0     0     0      0      0      0
TOTAL # 1kgSNPs   210   145   129   110    91    86    82    72    62     56     51     44
bio.13 bio.14 bio.15 bio.16 bio.17 bio.18 bio.19 bio.20 bio.21 bio.22 bio.23
rs10088218           9      8      8      8      8      6      1      1      0      0      0
rs11651755           0      0      0      0      0      0      0      0      0      0      0
rs146746174          2      1      0      0      0      0      0      0      0      0      0
rs17631303           3      1      0      0      0      0      0      0      0      0      0
rs1802669            1      1      1      1      1      1      1      1      1      1      1
rs4631563            0      0      0      0      0      0      0      0      0      0      0
rs4691139            1      1      1      1      1      1      1      1      1      1      1
rs4808075            1      1      1      1      1      0      0      0      0      0      0
rs61776206           0      0      0      0      0      0      0      0      0      0      0
rs62274042          19     16     13     10      9      5      3      1      1      1      1
rs6755777            0      0      0      0      0      0      0      0      0      0      0
rs7207826            1      1      0      0      0      0      0      0      0      0      0
rs76837345           0      0      0      0      0      0      0      0      0      0      0
TOTAL # 1kgSNPs     37     30     24     21     20     13      6      4      3      3      3
bio.24 bio.25
rs10088218           0      0
rs11651755           0      0
rs146746174          0      0
rs17631303           0      0
rs1802669            0      0
rs4631563            0      0
rs4691139            0      0
rs4808075            0      0
rs61776206           0      0
rs62274042           1      1
rs6755777            0      0
rs7207826            0      0
rs76837345           0      0
TOTAL # 1kgSNPs      1      1"



# Transformando todas essas informações em anotações de biofeatures, tudo 
# em formato .bed, essas são as anotações que integram os dados dos SNPs
# nas caracteristicas biologicas de cada cromossomo - biofeature
library(FunciSNP)

# Essa vai anotar todo YAFSNP identificado por ele de distância conhecida mais próxima de
# todas as características genômicas conhecidas (exon, intron, 5'UTR, 3'UTR,
# promotor, lincRNA ou no gene deserto (intergentic)) são usados  
# para anotar cada YAFSNP recentemente identificados como descrito acima.
# argumentos: um snp.list ==> ovarian. LEMBRANDO QUE: um objeto snp.list representa a saida de FunciSNP da função getFSNP()
# assim armazena em ovarian.anno todas as informações geradas em ovarian da outra função getFSNP(), em formato de data.frame, que será usada em outras funções
# Isso vai anotar tudo YAFSNP identificado por ele de distância da conhecida mais próxima TSS, se overlapps um 
# exon conhecida, intron, 5'UTR, 3'UTR, promotor, lincRNA ou no gene do deserto (intergentic) regiões.
ovarian.anno <- FunciSNPAnnotateSummary(ovarian)
class(ovarian.anno)
"[1] data.frame"

# transformando todas essas informações em anotações, tudo em formato .bed
# essas são as anotações que integram os dados dos SNPs nas caracteristicas
# biologicas de cada cromossomo - biofeature
ovarian.anno
"
  chromosome bio.feature.start
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep1_1.filtered          5           1293255
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep2_1.filtered          5           1293296
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep1_1.filtered             5           1293385
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep2_1.filtered             5           1293293
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep1_1.filtered              5           1293361
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep2_1.filtered              5           1293300
bio.feature.end
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep1_1.filtered         1295292
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep2_1.filtered         1295242
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep1_1.filtered            1295322
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep2_1.filtered            1295301
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep1_1.filtered             1295254
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep2_1.filtered             1295264
bio.feature
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep1_1.filtered iFTSEC246_H3K27ac_NA_rep1_1.filtered
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep2_1.filtered iFTSEC246_H3K27ac_NA_rep2_1.filtered
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep1_1.filtered       iOSE11_H3K27ac_NA_rep1_1.filtered
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep2_1.filtered       iOSE11_H3K27ac_NA_rep2_1.filtered
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep1_1.filtered         iOSE4_H3K27ac_NA_rep1_1.filtered
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep2_1.filtered         iOSE4_H3K27ac_NA_rep2_1.filtered
corr.snp.id
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep1_1.filtered chr5:1294263
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep2_1.filtered chr5:1294263
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep1_1.filtered    chr5:1294263
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep2_1.filtered    chr5:1294263
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep1_1.filtered     chr5:1294263
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep2_1.filtered     chr5:1294263
corr.snp.position tag.snp.id
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep1_1.filtered           1294263 rs10069690
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep2_1.filtered           1294263 rs10069690
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep1_1.filtered              1294263 rs10069690
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep2_1.filtered              1294263 rs10069690
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep1_1.filtered               1294263 rs10069690
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep2_1.filtered               1294263 rs10069690
tag.snp.position D.prime
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep1_1.filtered          1279790      NA
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep2_1.filtered          1279790      NA
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep1_1.filtered             1279790      NA
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep2_1.filtered             1279790      NA
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep1_1.filtered              1279790      NA
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep2_1.filtered              1279790      NA
R.squared p.value
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep1_1.filtered        NA       1
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep2_1.filtered        NA       1
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep1_1.filtered           NA       1
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep2_1.filtered           NA       1
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep1_1.filtered            NA       1
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep2_1.filtered            NA       1
distance.from.tag
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep1_1.filtered             14473
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep2_1.filtered             14473
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep1_1.filtered                14473
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep2_1.filtered                14473
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep1_1.filtered                 14473
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep2_1.filtered                 14473
population.count population
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep1_1.filtered              379        EUR
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep2_1.filtered              379        EUR
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep1_1.filtered                 379        EUR
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep2_1.filtered                 379        EUR
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep1_1.filtered                  379        EUR
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep2_1.filtered                  379        EUR
nearest.lincRNA.ID
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep1_1.filtered     TCONS_00010241
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep2_1.filtered     TCONS_00010241
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep1_1.filtered        TCONS_00010241
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep2_1.filtered        TCONS_00010241
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep1_1.filtered         TCONS_00010241
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep2_1.filtered         TCONS_00010241
nearest.lincRNA.distancetoFeature
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep1_1.filtered                           -132845
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep2_1.filtered                           -132845
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep1_1.filtered                              -132845
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep2_1.filtered                              -132845
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep1_1.filtered                               -132845
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep2_1.filtered                               -132845
nearest.lincRNA.coverage
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep1_1.filtered                 upstream
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep2_1.filtered                 upstream
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep1_1.filtered                    upstream
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep2_1.filtered                    upstream
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep1_1.filtered                     upstream
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep2_1.filtered                     upstream
nearest.TSS.refseq
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep1_1.filtered          NM_198253
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep2_1.filtered          NM_198253
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep1_1.filtered             NM_198253
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep2_1.filtered             NM_198253
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep1_1.filtered              NM_198253
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep2_1.filtered              NM_198253
nearest.TSS.GeneSymbol
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep1_1.filtered                   TERT
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep2_1.filtered                   TERT
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep1_1.filtered                      TERT
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep2_1.filtered                      TERT
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep1_1.filtered                       TERT
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep2_1.filtered                       TERT
nearest.TSS.ensembl
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep1_1.filtered     ENST00000310581
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep2_1.filtered     ENST00000310581
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep1_1.filtered        ENST00000310581
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep2_1.filtered        ENST00000310581
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep1_1.filtered         ENST00000310581
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep2_1.filtered         ENST00000310581
nearest.TSS.coverage
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep1_1.filtered               inside
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep2_1.filtered               inside
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep1_1.filtered                  inside
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep2_1.filtered                  inside
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep1_1.filtered                   inside
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep2_1.filtered                   inside
nearest.TSS.distancetoFeature
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep1_1.filtered                           899
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep2_1.filtered                           899
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep1_1.filtered                              899
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep2_1.filtered                              899
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep1_1.filtered                               899
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep2_1.filtered                               899
Promoter utr5 Exon Intron
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep1_1.filtered       NO   NO  YES    YES
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep2_1.filtered       NO   NO  YES    YES
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep1_1.filtered          NO   NO  YES    YES
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep2_1.filtered          NO   NO  YES    YES
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep1_1.filtered           NO   NO  YES    YES
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep2_1.filtered           NO   NO  YES    YES
utr3 Intergenic
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep1_1.filtered   NO         NO
rs10069690:EUR.chr5:1294263.iFTSEC246_H3K27ac_NA_rep2_1.filtered   NO         NO
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep1_1.filtered      NO         NO
rs10069690:EUR.chr5:1294263.iOSE11_H3K27ac_NA_rep2_1.filtered      NO         NO
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep1_1.filtered       NO         NO
rs10069690:EUR.chr5:1294263.iOSE4_H3K27ac_NA_rep2_1.filtered       NO         NO
"
# esse corte faz parte da "chamada de picos" - enriquecimento dos dados de
# sequenciamento é o ponto de corte, ou seja começa o corte a partir deste
# numero, digitado na função é daí em diante a busca para mostrar na tabela..
# ==> que ainda está na regiao de chamada de picos..
# Mostra a quantidade total na região, o corte, e a porcentagem: de tagSNPs,
# de SNPs e de biofeatures.[ie, a integração]
# BASICAMENTE:  Utilizando um valor especificado Rsquare (0-1) para o subconjunto dos 
# dados, é gerado uma tabela que RESUME O NÚMERO TOTAL DE YAFSNPS, TAGSNPS ASSOCIADO, E NÚMERO DE BIOFEATURES SOBREPOSTAS.
#  ISTO IRÁ FORNECER USUÁRIO UM PRIMEIRO OLHAR PARA O NÚMERO TOTAL DE YAFSNP DISPONÍVEL NO PONTO DE CORTE RSQUARE PARTICULAR.
FunciSNPtable(ovarian.anno, rsq=0.8)

# print ovarian
"
Total R.sq>=0.8 Percent
tagSNPs        19        16   84.21
1K SNPs     12532       208    1.66
Biofeatures    35        32   91.43
"
# especifica os pontos de cortes dos genes mais específicos: assim é 
# mostrado uma tabela com os nomes dos genes em que esse rsq corta
# => INFORMA OS GENES MAIS PRÓXIMOS DOS YAFSNPS DESSES PONTOS DE CORTES
FunciSNPtable(ovarian.anno, rsq=0.8, geneSum=TRUE)
"
 Gene_Names
1       ABHD8
2         ABO
3      ANKLE1
4    ARHGAP27
5        BNC2
6     C4orf39
7      CHMP4C
8      DNALI1
9        GNL2
10       GPX6
11      HNF1B
12      HOXD3
13      HOXD4
14  LOC401022
15  LOC730091
16    MIR1208
17     MLLT10
18    PA2G4P4
19    PLEKHM1
20      SKAP1
21      SNIP1
22     TIPARP
23 TIPARP-AS1
"

# RESUMO DA SOBREPOSIÇÃO DE SNPS EM BIOFEATURES
# nesta função ocorre um resumo das sobreposições recentemente
# descobertas -> YAFSNPs
FunciSNPsummaryOverlaps(ovarian.anno)  # argumento: data.frame

# print dos dados do ovarian.anno:
"

                bio.1 bio.2 bio.3 bio.4 bio.5 bio.6 bio.7 bio.8 bio.9 bio.10 bio.11 bio.12
rs10069690        248    99    69    59    52    38    29    21    16     11      8      4
rs10088218        535   174   134   109    92    81    76    66    65     54     52     40
rs11651755        478   346   274   227   194   168   155   130   109     98     85     70
rs146746174       185   129   105    96    79    76    73    69    59     54     49     45
rs17631303        225   156   121   103    86    72    68    62    49     43     36     31
rs1802669         149    98    78    70    62    61    56    54    49     43     35     31
rs2165805         518   234   185   148   127   109    97    84    81     67     63     49
rs2191035         175   109    87    76    61    56    51    44    37     34     26     24
rs2519093         319   181   151   132   108    83    71    60    49     43     32     29
rs4631563         250   188   136   109    83    73    65    64    64     53     47     42
rs4691139         196    93    66    54    48    43    39    39    38     30     26     17
rs4808075         478   377   324   287   266   237   215   193   177    155    136    126
rs61776206        273   198   159   138   116   106    97    95    90     86     79     70
rs62274042        561   506   463   410   364   337   306   266   246    223    195    178
rs6755777         360   314   269   223   181   156   122    93    70     56     50     42
rs7207826         308   238   201   169   143   126   105    96    86     75     66     60
rs757210          483   351   278   229   195   169   155   130   109     98     85     70
rs7671665         444   360   286   246   203   174   160   140   125    109     93     75
rs76837345        221   133   113    93    83    66    54    47    44     39     26     19
TOTAL # 1kgSNPs  6406  4284  3499  2978  2543  2231  1994  1753  1563   1371   1189   1022

               
                bio.13 bio.14 bio.15 bio.16 bio.17 bio.18 bio.19 bio.20 bio.21 bio.22 bio.23
rs10069690           4      4      3      2      0      0      0      0      0      0      0
rs10088218          34     30     28     28     27     20     13     11      7      5      3
rs11651755          51     42     37     32     29     20     19     10      7      5      2
rs146746174         38     26     21     17     15     11     10      7      0      0      0
rs17631303          27     21     14     11     10      9      5      3      2      2      2
rs1802669           29     27     22     19     17     15      9      7      5      3      2
rs2165805           43     38     35     34     32     25     17     15      9      5      3
rs2191035           21     17     11      6      5      2      2      2      2      2      0
rs2519093           24     20     15      4      2      1      0      0      0      0      0
rs4631563           37     30     23     14     10      5      3      3      2      2      0
rs4691139           13     13     12     11     10      8      4      3      2      2      2
rs4808075          107    102     87     71     47     32     16      3      2      0      0
rs61776206          61     52     40     31     22     14      6      4      2      0      0
rs62274042         159    140    120    100     90     68     52     33     13      6      4
rs6755777           33     21     15      9      1      0      0      0      0      0      0
rs7207826           55     53     48     39     32     25     16     13     12      8      8
rs757210            51     42     37     32     29     20     19     10      7      5      2
rs7671665           61     56     40     35     32     24     15     11     10      7      4
rs76837345          11     10      9      8      7      5      2      0      0      0      0
TOTAL # 1kgSNPs    859    744    617    503    417    304    208    135     82     52     32


                bio.24 bio.25 bio.26 bio.27 bio.28 bio.29
rs10069690           0      0      0      0      0      0
rs10088218           3      3      2      2      2      1
rs11651755           2      0      0      0      0      0
rs146746174          0      0      0      0      0      0
rs17631303           2      2      1      1      0      0
rs1802669            0      0      0      0      0      0
rs2165805            3      3      2      2      2      1
rs2191035            0      0      0      0      0      0
rs2519093            0      0      0      0      0      0
rs4631563            0      0      0      0      0      0
rs4691139            0      0      0      0      0      0
rs4808075            0      0      0      0      0      0
rs61776206           0      0      0      0      0      0
rs62274042           3      2      0      0      0      0
rs6755777            0      0      0      0      0      0
rs7207826            5      3      1      0      0      0
rs757210             2      0      0      0      0      0
rs7671665            4      3      2      1      0      0
rs76837345           0      0      0      0      0      0
TOTAL # 1kgSNPs     24     16      8      6      4      2

"

# RESUMO DA SOBREPOSIÇÃO DE SNPS EM BIOFEATURES, A PARTIR DO CORTE ESPECIFICADO
# Essa ultima função mostra um resumo dos YAFSNPs nos cromossomos, a partir
# dos corte especificado, na função --> rsq=0.44
# MESMA FUNÇÃO DE CIMA, MAS MOSTRA NUM CORTE ESPECIFICADO
FunciSNPsummaryOverlaps(ovarian.anno, rsq=0.8)
"

                bio.1 bio.2 bio.3 bio.4 bio.5 bio.6 bio.7 bio.8 bio.9 bio.10 bio.11 bio.12
rs10088218         44    27    22    19    18    18    18    18    18     18     17     13
rs11651755          4     4     3     2     2     2     2     2     0      0      0      0
rs146746174        46    26    21    18    13    12    12    10     7      5      5      4
rs17631303         50    28    23    19    13    12    12    10     7      5      5      4
rs1802669          10     7     6     6     4     4     4     4     3      2      1      1
rs2191035           1     0     0     0     0     0     0     0     0      0      0      0
rs2519093           2     0     0     0     0     0     0     0     0      0      0      0
rs4631563           2     1     0     0     0     0     0     0     0      0      0      0
rs4691139           1     1     1     1     1     1     1     1     1      1      1      1
rs4808075           5     4     4     3     2     2     2     2     1      1      1      1
rs61776206          6     3     1     1     1     1     1     1     1      1      1      1
rs62274042         60    57    55    48    42    39    36    30    30     29     26     23
rs6755777          12    10    10     7     4     3     0     0     0      0      0      0
rs7207826           7     5     3     3     1     1     1     1     1      1      1      1
rs757210            2     0     0     0     0     0     0     0     0      0      0      0
rs76837345          2     1     1     0     0     0     0     0     0      0      0      0
TOTAL # 1kgSNPs   254   174   150   127   101    95    89    79    69     63     58     49


                bio.13 bio.14 bio.15 bio.16 bio.17 bio.18 bio.19 bio.20 bio.21 bio.22 bio.23
rs10088218          10      9      9      9      9      6      1      1      0      0      0
rs11651755           0      0      0      0      0      0      0      0      0      0      0
rs146746174          3      1      0      0      0      0      0      0      0      0      0
rs17631303           3      1      0      0      0      0      0      0      0      0      0
rs1802669            1      1      1      1      1      1      1      1      1      1      1
rs2191035            0      0      0      0      0      0      0      0      0      0      0
rs2519093            0      0      0      0      0      0      0      0      0      0      0
rs4631563            0      0      0      0      0      0      0      0      0      0      0
rs4691139            1      1      1      1      1      1      1      1      1      1      1
rs4808075            1      1      1      1      1      0      0      0      0      0      0
rs61776206           1      1      0      0      0      0      0      0      0      0      0
rs62274042          20     17     14     11     10      6      4      1      1      1      1
rs6755777            0      0      0      0      0      0      0      0      0      0      0
rs7207826            1      1      0      0      0      0      0      0      0      0      0
rs757210             0      0      0      0      0      0      0      0      0      0      0
rs76837345           0      0      0      0      0      0      0      0      0      0      0
TOTAL # 1kgSNPs     41     33     26     23     22     14      7      4      3      3      3


                bio.24 bio.25
rs10088218           0      0
rs11651755           0      0
rs146746174          0      0
rs17631303           0      0
rs1802669            0      0
rs2191035            0      0
rs2519093            0      0
rs4631563            0      0
rs4691139            0      0
rs4808075            0      0
rs61776206           0      0
rs62274042           1      1
rs6755777            0      0
rs7207826            0      0
rs757210             0      0
rs76837345           0      0
TOTAL # 1kgSNPs      1      1
"

# DIGITANDO A MESMA FUNÇÃO ANTERIORMENTE PARA CARREGAR NOVAMENTE E NAO DAR ERRO!!!
ovarian.anno <- FunciSNPAnnotateSummary(ovarian)

# MOSTRA OS GRAFICOS!!! VARIOS TIPOS, PARA ESCOLHER A MELHOR MANEIRA DE SE REPRESENTAR OS DADOS GERADOS!!!!!!
# FunciSNPplot:: PARA VISUALIZAR RESUMO DOS YAFSNP.
# mostra a quantidade de YAFSNPs em cada ponto de corte, i e, rsquare 
# relacionado --> em cada ponto contem um numero total de YAFSNPs
# correlacionado
FunciSNPplot(ovarian.anno)   # MOSTRAR FIGURA: plot1; argumento: data.frame

# mostra o TOTAL DE YAFSNPS DIVIDIDO PELO TAGSNP ASSOCIADO, DIVIDE A DISTRIBUIÇÃO DE SNPs POR TAGSNP.
#(É MAIS ESPECIFICADO)
FunciSNPplot(ovarian.anno, splitbysnp=TRUE)  # MOSTRAR FIGURA: plot2; argumento: data.frame

# salvando o arquivo em pdf - a imagem - plot
ggsave("ovarian_dist_bysnp.pdf")

# HEATMAP:: CORRELAÇÃO PARA VISUALIZAR O NÚMERO DE SNPS CORRELACIONADOS EM CADA TAGSNP SOBREPONDO 
# CADA CARACTERÍSTICA BIOLÓGICA. MAIS INFORMATIVO SE USADO COM UM VALOR RSQ.
# mostra a relação de biofatures e tagSNPs associados a uma quantidade de 
# SNPs, neste corte.
FunciSNPplot(ovarian.anno, heatmap=TRUE, rsq=0.8, heatmap.key=TRUE)  # MOSTRAR FIGURA: plot3; argumento: data.frame

# Este tipo de plot informa o ENRIQUECIMENTO RELATIVO A CARACTERÍSTICAS GENÔMICAS
# genomicSum empilhadas GRÁFICO DE BARRAS RESUMINDO TODOS OS SNPS CORRELACIONADOS PARA CADA UMA DAS CARACTERÍSTICAS 
# GENOMICAS IDENTIFICADOS (EXON, INTRON, 5'UTR, 3'UTR, PROMOTOR, LINCRNA OU NO GENE DESERTO (INTERGENTIC)). MAIS INFORMATIVO SE USADO COM UM VALOR RSQ.
# com todos os SNPs e do outro lado num corte maior ou igual a 
# 0.5 de rsquare.
FunciSNPplot(ovarian.anno, rsq=0.8, genomicSum=TRUE, save=F)  # MOSTRAR FIGURA: plot4; argumento: data.frame

# essa função ajuda a IDENTIFICAR DE UMA MANEIRA MELHOR, OS TAGSNPS 
# ASSOCIADOS AOS BIOFEATURES E TAMBÉM OS YAFSNPS, DE UMA MANEIRA MAIS 
# ILUSTRATIVA. SALVANDO EM FORMATO ".bed". 
FunciSNPbed(ovarian.anno, rsq=0.8)  # MOSTRAR O ARQUIVO EM FORMATO ".bed" E NO SITE BROWSER
"
####
Bed file "FunciSNP_results_rsq.0.8.bed" created successfully.
(See folder: "/home/mikely")
Total corSNP (RED):  208 
Total tagSNP (BLK):  16 

To view results, submit bed file as a
custom track in UCSC Genome Browser (genome.ucsc.edu), 

Now have fun with your new YAFSNPss!!
####
"













