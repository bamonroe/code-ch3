
ybot <- 0
ytop <- 90

xvline <- .1





p <- ggplot(data = data.frame(x = 0), mapping = aes(x=x)) +
    layer(stat = "function",
          fun = UA,
          mapping = aes(color = "UA")
          ) +
    layer(stat = "function",
          fun = UBuy,
          mapping = aes(color = "UBuy")
          ) +
    scale_x_continuous(
    				   limits = c(-.25, 1),
    				   breaks=c(seq(from=-.2, to=1, by=.15))
    				   ) +
    scale_y_continuous(
    				   limits = c(ybot, ytop),
    				   breaks=c(seq(from=ybot, to=ytop, by=15))
    				   ) +
    scale_color_manual(
    				   name = "Functions",
                       values = c("blue", "red"), # Color specification
                       labels = c("Utility of Ticket", "Utility of Buy Price","Utility of Sell Price")
                       ) +
	theme(
		axis.text.x=element_text(angle=50, size=12, vjust=0.5), 
		axis.title.x = element_text(color="forestgreen", vjust=-0.35),
		axis.title.y = element_text(color="forestgreen" , vjust=0.35)
		) +
	labs(x="CRRA Value", y="Utility", title="Buying the Ticket")+ 
	geom_segment(aes(x = xvline, y = ybot, xend = xvline, yend = UA(xvline)), , linetype="longdash") 
	geom_segment(aes(x = xvline, y = UA(xvline), xend = xvline, yend = UA(xvline)), , linetype="longdash") +




p


q <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) +
   layer(stat = "function",
         fun = CEA,
         mapping = aes(color = "UA")
         ) +
   layer(stat = "function",
         fun = CEBuy,
         mapping = aes(color = "UBuy")
         ) +
   layer(stat = "function",
         fun = CEC,
         mapping = aes(color = "USell")
         ) +
   scale_x_continuous(
   				   limits = c(-.25, 1)
   				   ) +
   scale_y_continuous(
   				   limits = c(25, 65)
   				   ) +
   scale_color_manual(
   				   name = "Functions",
                      values = c("blue", "red","black"), # Color specification
                      labels = c("Utility of Ticket", "Utility of Buy Price","Utility of C")
                      ) +
   theme(
   	axis.text.x=element_text(angle=50, size=12, vjust=0.5), 
   	axis.title.x = element_text(color="forestgreen", vjust=-0.35),
   	axis.title.y = element_text(color="forestgreen" , vjust=0.35)
   	) +
   labs(x="CRRA Value", y="Utility", title="Buying the Ticket")+
   geom_vline(xintercept=0, linetype="longdash")
q


cat("\014")
cat("\014")
