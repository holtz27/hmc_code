library(ggplot2)
library(cowplot)
library(zoo)
library(moments)

#Banco coletado em: 
#https://br.investing.com/indices/bovespa-historical-data
#Dados de 04.01.2010 a 30.12.2019

ibovespa = read.csv('Ibovespa Dados Históricos.csv')
ibovespa$Data = as.Date(ibovespa$Data, format = '%d.%m.%Y')
#Selecionando subconjunto
data.ini = '2010-01-04'
data.fim = '2019-12-30'
ibovespa = subset(ibovespa, Data >=  data.ini & Data <= data.fim)
#View(ibovespa)

#indices ibovespa
indice = ibovespa$Último
indice = indice[length(indice):1]

#datas
dia = ibovespa$Data
dia = dia[length(dia):1]

#log - retornos r = log(pt/pt-1) calculados
ret = indice[2:length(indice)]/indice[-length(indice)]
log.ret = log(ret)

#Gerando data frames
df.1 = data.frame(dia, indice)
df.2 = data.frame(x=dia[-1], y=log.ret)

#Gerando figuras
#Figura preços
g = ggplot(df.1, aes(dia,indice)) + geom_line()
g = g + ylab('Ìndice') + xlab('Dias') + theme_bw(base_size = 18)
g = g + scale_x_date(date_breaks = "2 year", date_labels = "%Y")
g = g + theme(plot.caption = element_text(size = 15, hjust = 0.5),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_text(size=14),
              axis.title.x = element_text(vjust = -2),
              axis.text = element_text(colour = "black"))
g

#Figura retornos
h = ggplot(df.2, aes(x,y)) + geom_line()
h = h + ylab('Retornos') + xlab('Dias') + theme_bw(base_size = 18)
h = h + scale_x_date(date_breaks = "2 year", date_labels = "%Y")
h = h + theme(plot.caption = element_text(size = 15, hjust = 0.5),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_text(size=14),
              axis.title.x = element_text(vjust = -2),
              axis.text = element_text(colour = "black"))
h

#Histograma retornos
i = ggplot(df.2, aes(x=y)) + geom_histogram(aes(y = ..density..), 
                                          fill="grey", 
                                          color="white")
#i = i + stat_function(fun = dnorm, args = list(mean = mean(log.ret),
#                                               sd = sd(log.ret)),
#                      linetype = "dashed")
i = i + ylab('Frequência') + xlab('') + theme_bw(base_size = 18)
i = i + theme(plot.caption = element_text(size = 15, hjust = 0.5),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_text(size=14),
              axis.title.x = element_text(vjust = -2),
              axis.text = element_text(colour = "black"))
#i = i + geom_density(adjust = 1.5)
i

#QQ plot retornos
j = ggplot(df.2, aes(sample = y)) + geom_qq()
j = j + stat_qq_line()
j = j + ylab('Quantis amostrais') + xlab('Quantis teóricos') + theme_bw(base_size = 18)
j = j + theme(plot.caption = element_text(size = 15, hjust = 0.5),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_text(size=14),
              axis.title.x = element_text(vjust = -2),
              axis.text = element_text(colour = "black"))
j

#Descritiva dos log-retornos
mean(log.ret)         #Média
median(log.ret)       #Mediana
sd(log.ret)           #Desvio padrão
skewness(log.ret)     #Assimetria
kurtosis(log.ret)     #Curtose
min(log.ret)          #Mínimo
max(log.ret)          #Máximo

