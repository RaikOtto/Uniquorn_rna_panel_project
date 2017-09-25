library(stringr)

t = read.table("~/Downloads/spia_cls.csv",sep="\t",header = T)

spia_cls = str_replace_all(t$Spia_cls, pattern = "-","")
spia_cls = str_replace_all(spia_cls, pattern = " ","")
spia_cls = str_replace_all(spia_cls, pattern = "\\.","")
spia_cls = str_to_upper( spia_cls  )
spia_cls = spia_cls[ spia_cls!= ""]

spia_cls[ spia_cls == "U87" ] = "U87MG"
spia_cls[ spia_cls == "U138" ] = "U138MG"
spia_cls[ spia_cls == "H128" ] = "NCIH128"
spia_cls[ spia_cls == "H1339" ] = "NCIH1339"
spia_cls[ spia_cls == "H1395" ] = "NCIH1395"
spia_cls[ spia_cls == "H1437" ] = "NCIH1437"
spia_cls[ spia_cls == "H1450" ] = "NCIH1450"
spia_cls[ spia_cls == "H1770" ] = "NCIH1770"
spia_cls[ spia_cls == "H1819" ] = "NCIH1819"
spia_cls[ spia_cls == "H2009" ] = "NCIH2009"
spia_cls[ spia_cls == "H2141" ] = "NCIH2141"
spia_cls[ spia_cls == "H2171" ] = "NCIH2171"
spia_cls[ spia_cls == "H2195" ] = "NCIH2195"
spia_cls[ spia_cls == "H220" ]  = "NCIH220"
spia_cls[ spia_cls == "H378" ]  = "NCIH378"
spia_cls[ spia_cls == "H460" ]  = "NCIH460"
spia_cls[ spia_cls == "H838" ]  = "NCIH838"
spia_cls[ spia_cls == "H889" ]  = "NCIH889"

index = which( spia_cls %in% t$Name_query)
spia_cls[index]
spia_cls[-index]

res_table = t
res_table$Spia_cls = spia_cls

write.table(spia_cls,"~/Downloads/res_table.tab",sep = "\t", quote = F, row.names = F)
t = read.table("~/Downloads/res_table.tab",sep="\t",header = T)

table( t$Identification_successful[ which( t$Name_query %in%  spia_cls ) ] )
