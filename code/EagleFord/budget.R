
require("ggplot2")
dat <- data.frame(money=c(1, 0.5, 0.5, 1.5), 
                     Budget=c("$1 mln Vendor contracts", 
                            "$0.5 mln software&compute charges",
                            "1 FTE PTD/TASE(Statistics&Chemometrics",
                            "3 FTEs PTI/RP, PTI/UU and PTI/C")
                     )

#pie(budget)


p <- ggplot(dat, aes(x=1, y=money, fill=Budget)) +
  geom_bar(stat="identity") +
  ggtitle(" ")

p <- p + coord_polar(theta='y')
print(p)




