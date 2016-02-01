
## Carga librerías
# librería para leer datos epa
library(MicroDatosEs)

library(plyr)

# librerías para leer shp, mapas, etcmapas 
library(maptools)
library(rgdal)
library(sp)
library(RColorBrewer)
library(classInt)
library(ggmap)
library(ggplot2)
library(scales)

# otras librerías
library(parallel)
library(lme4)
library(boot)

## bajamos los datos de mi dropbox, podemos cambiar method por "wget" "auto" 
# si no funciona con "curl", ver ayuda de download.file

download.file("https://dl.dropboxusercontent.com/u/2712908/EPA/EPAT0415",
              destfile="EPAT0415", method="curl")

# epa <- epa2005("EPA/EPAT0415")
epa <- epa2005("EPAT0415")


# Seleccionamos solo las variables que nos interesan,provincia, edad, estudios, aoi
# que indica la situación laboral, (activo, inactivo, parado, ocupado) y el factor
# de elevación ya que la epa está ponderada

dat <- subset(epa, select=c(prov,edad,nforma, aoi,factorel))

# borramos objeto epa
rm(epa) 


# recodificamos usando recode de la función memisc, ya que dat es es un objeto de tipo 
# data.set del paquete memisc class(dat)

recodificacion <- function (dat) {
   dat$aoi <- memisc::recode(dat$aoi, "o" = 1 <- 3:4, "p" = 2 <- 5:6, "i" = 3 <- 7:9)
   
   dat$nforma3 <- memisc::recode(dat$nforma,                        
                                 "Est primarios o menos"  = 1 <- c("AN","P1","P2"),
                                 "Est. Secundarios" = 2 <- c("S1","SG","SP"),
                                 "Est. Universitarios "  = 3 <- c("SU")
   )
   dat$gedad <- memisc::recode(dat$edad,
                               "15 años o menos " = 1 <- c(0,5,10),
                               "De 16 a 34 " = 2 <- c(16,20,25,30),
                               "De 35 a 54" = 3 <- c(35,40,45,50),
                               "De 55 o más" = 4 <- c(55,60,65)
   )
   
   dat
}


dat <- recodificacion(dat)

# convierto a data.frame de toda la vida

dat <- as.data.frame(dat)

# variable peso 
dat$peso <- dat$factorel*nrow(dat)/sum(dat$factorel)


# creo variable numérica para las provincias
dat$prov.n <- as.numeric(dat$prov)

# eliminar menores de 16 años  
dat <- dat[ as.numeric(dat$edad) > 3, ]

# eliminar inactivos
dat <- dat[ dat$aoi != "i", ]

# eliminar niveles que no se usan
dat$gedad <- droplevels( dat$gedad)
levels(dat$gedad)


# Función que ajusta un modelo mixto sencillo y que calcule
# las medias de las probabilidades ajustadas en cada combinación de prov, gedad y nforma3
predicciones <- function( data, indices) {
   d <- data[indices, ]
   fit <- glmer(aoi=="p" ~ (1|prov.n) + (1|gedad) + (1|nforma3), family=binomial
                ,data=d, nAGQ = 0)
   
   pred <- fitted(fit)
   return(tapply(pred, list(d$prov.n,d$gedad, d$nforma3),mean, na.rm=TRUE))
}


# Bootstrap, suele tardar bastante (500 bootstrap más de 3 horas)
# pendiente de probar utilizar brms, que usa stan y poner  a prioris informativas
# sacadas de la estimación de los efectos aleatorios en anteriores epas.
# he comprobado que si en vez de 500 hacemos 100 o menos
# las diferencias del ajuste son mínimas.

# en el bootstrap ajustamos weights = peso para tener en cuenta la ponderación de la encuesta.
# y que cada muestra bootstrap sea representativa de la población. (toma muestras
# con reemplazamiento y con probabilidad de selección de cada unidad como 1/peso)

# en puridad habría que hacer el bootstrap dentro de cada provincia.

# Podría haber utilizado la función bootMer del paquete lme4, pero no recoge todas las 
# fuentes de variabilidad, así que he preferido remuestrear la epa y en cada muestra
# ajustar el modelo de nuevo. 

#--------------------------------------------------
# OJO, tarda mucho !! probar con menos n.bootstrap
#---------------------------------------------------

n.bootstrap <- 500
   system.time (
      res.boot <- boot(data = dat, statistic = predicciones, R = n.bootstrap, weights = dat$peso, parallel="multicore",
                    ncpus= detectCores()-1)
               )

head(res.boot$t)

boot.ci(res.boot, type = "perc")

# intervalos de confianza al 95% 

estim <- data.frame(fit = apply(res.boot$t, 2, function(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE))),
           lwr = apply(res.boot$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE))),
           upr = apply(res.boot$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE))))


# de 1 a 52 es primarios o menos de 16 a 34 según cada provincia, 
# de 53 a 104 son primarios o menos de 35 a 54 en cada provincia
# de 105 a 156 son primarios o menos de más de 54 
# en cada provincia



estudios <- c(rep("Est primarios o menos",52*3), rep("Est. Secundarios",52*3),rep("Est Universitarios",52*3))
edades <- rep(rep(c("De 16 a 34","De 35 a 54","De 55 o más"),each=52),3)

provincias <- as.data.frame(table(dat$prov.n,dat$prov))
provincias <- provincias[provincias$Freq>0,]

prov.nombres <- rep(provincias$Var2,9)
prov.n <-  rep(provincias$Var1,9)



estim$provincia <- prov.nombres
estim$estudios <- estudios
estim$gedad <- edades
estim$prov <- prov.n

# prov en el data.frame estim es numérica, para unir bien en el mapa tenemos que añadir un 0 a
# los valores con un solo dígito

estim$prov <- as.character(estim$prov)

estim$prov <- ifelse(nchar(estim$prov)<2, paste0(0,estim$prov), estim$prov)

# save(estim,file = "estimaciones4t2005.RData")

# ejemplo para estudios universitarios con intervalos de confianza
ggplot(estim[estim$estudios=="Est Universitarios" , ], aes(x=reorder(provincia,fit), y=fit))+
   geom_point() + 
   geom_errorbar(aes(ymax=upr, ymin=lwr, width= 0.2))+
   facet_wrap( ~gedad) +
   scale_y_continuous(labels=percent) +
   theme(axis.text.y = element_text(size=rel(1.4)),
      axis.text.x = element_text(size = rel(1.4)),
      plot.title = element_text(face="bold", size=rel(1.4)),
      legend.text = element_text(size=rel(1.4)),
      strip.text = element_text(size=rel(1)),
      axis.line = element_line(colour = "black")
   
      # ,
      # panel.grid.major = element_blank()
      # ,
      # panel.grid.minor = element_blank()
      # ,
      # panel.border = element_blank(),
      # panel.background = element_rect(fill="grey97")
   ) +
      
   labs(list(x="",y="")) +
   coord_flip()+
   ggtitle("Tasa de paro \npor edad para población\n con estudios universitarios")



##########################################
# MAPA todas provincias
##########################################


# bajamos la capa con las provincias, (tb está en la página del INE)

download.file("https://dl.dropboxusercontent.com/u/2712908/spain_provinces_ind_2.zip",
              destfile = "spain_provinces_ind2.zip", method="auto")

unzip("spain_provinces_ind2.zip")


# leemos con readOGR que convierte bien los strings

prov.map <- readOGR(dsn=getwd(),layer="spain_provinces_ind_2")

## ------------------------------------------------------------------------
# En el objeto prov.map tenemos una variable con el nombre de las provincias
prov.map$NOMBRE99
# pero es mejor quedarse con variable del código
prov.map$PROV
# pongo en minúsculas las variables en prov.map@data
colnames(prov.map@data) <- tolower(colnames(prov.map@data))

## pasamos la capa a un objeto que puede tratar ggplot2 usando fortify e 
# indicaando la variable de unión con los datos 

prov.ggmap <- fortify(prov.map, region="prov")

# Hacemos intervalos

intervalos <- classIntervals(estim$fit,6,style="pretty")
intervalos

estim$intervalos <- as.factor(findInterval(estim$fit, intervalos$brks,all.inside=TRUE))

ggplot(estim) +
   geom_map(aes(map_id = prov, fill=intervalos),
            map = prov.ggmap,
            colour="grey30") +
   expand_limits(x = prov.ggmap$long, y = prov.ggmap$lat) + 
   facet_grid(  gedad ~ estudios)  + 
   scale_fill_brewer(palette="Reds",
                     labels=c("Menos del 10%","[10% - 20%)","[20% - 30%)",
                              "[30% - 40%)", "[40% - 50%)", "[50% - 60%)",
                              "[60% - 70%)"
                     )) +
   
   scale_x_continuous(breaks=NULL) +
   scale_y_continuous(breaks=NULL) +
   theme(axis.text.y = element_blank(),
         axis.text.x = element_blank(),
         plot.title = element_text(face="bold", size=rel(1.5)),
         legend.text = element_text(size=rel(1.5)),
         panel.background=element_rect(fill="#EEEDF6"),
         strip.text = element_text(face="bold",size=rel(1.5)))+
   labs(list(x="",y="",
             fill="")) +
   ggtitle("Tasa de paro 4T 2015\npor edad y estudios\n con modelo mixto")


