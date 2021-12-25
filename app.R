require(shiny)
require(shinydashboard)
require(shinyjs)
require(e1071)
require(tidyverse)

# Put all functions here

#一.欧式看涨期权空头Delta中性对冲
#输入：
#X:执行价格
#S0：股票现价
#r:无风险利率（年）
#sig:波动率
#Ti：到期时间（周）
#dt:调整间隔（周）
#N:交易数量
#输出：
#S：股价模拟值
#TABLE：对冲表格
#CostNow：对冲成本现值
DeltaEuroCall <- function(X=50,S0=49,r=0.05,sig=0.2,
                          Ti=20,dt=1,N=1e5){
    dt <- max(1/70,dt)
    dt <- min(Ti,dt)
    #生成几何布朗运动随机数
    S <- c(S0,rep(NA,floor(Ti/dt)))
    for(j in 2:length(S)){
        S[j] <- S[j-1]*exp((r-0.5*sig^2)*dt/52+sig*rnorm(1,0,sqrt(dt/52)))
    }
    #提取调整节点
    Week <- seq(0,Ti,by=dt)
    Stock_Price <- S
    Delta <-(log(Stock_Price/X)+(r+0.5*sig^2)*(Ti-Week)/52)/(sig*sqrt((Ti-Week)/52))
    Delta <- pnorm(Delta)
    #计算第0周
    Sh <- Shares_Purchased <- Delta[1]*N
    Cost <- Shares_Purchased*Stock_Price[1]
    Cu <- Cumulative_Cash <- Cost
    In <- Interest <- Cumulative_Cash*(exp(r*(dt/52))-1)
    #计算剩余量
    for(j in 2:length(Delta)){
        Sh <- (Delta[j]-Delta[j-1])*N
        Co <- Sh*Stock_Price[j]
        Cu <- Cu+In+Co
        In <- Cu*(exp(r*(dt/52))-1)
        Shares_Purchased <- c(Shares_Purchased,Sh)
        Cost <- c(Cost,Co)
        Cumulative_Cash <- c(Cumulative_Cash,Cu)
        Interest <- c(Interest,In)
    }
    Interest[length(Delta)] <- NA
    #整理
    TABLE <- data.frame(
        Week,Stock_Price,Delta,Shares_Purchased,Cost,
        Cumulative_Cash,Interest
    )
    #对冲成本现值
    CostNow <- ifelse(Stock_Price[length(Delta)]>X,
                      round((Cu-N*X)*exp(-r*Ti/52),2),
                      round(Cu*exp(-r*Ti/52),2))
    return(list(S,TABLE,CostNow))
}
#二、欧式看涨期权的Black Scholes公式
Black_Scholes_Call <- function(X=50,S0=49,r=0.05,sig=0.2,
                              Ti=20,N=1e5){
    d1 <- (log(S0/X)+(r+0.5*sig^2)*Ti/52)/(sig*sqrt(Ti/52))
    d2 <- d1-sig*sqrt(Ti/52)
    Price <- -X*exp(-r*Ti/52)*pnorm(d2)+S0*pnorm(d1)
    Price <- round(Price*N,2)
    return(Price)
}
#三、利用二叉树对美式看跌期权定价
#输入：
#n:期数
#输出：
#1.Price：期权价格
#2.TABLE：二叉树表
BiAmerCall <- function(X=50,S0=50,r=0.1,sig=0.4,
                       Ti=0.4167,n=5,N=1e5){
    u <- exp(sig*sqrt(Ti/n))
    d <- exp(-sig*sqrt(Ti/n))
    p <- (exp(r*Ti/n)-d)/(u-d)
    edt <- exp(-r*Ti/n)
    #构造数据结构
    result1 <- matrix(rep(NA,(n+1)^2),ncol=n+1)
    result2 <- matrix(rep(NA,(n+1)^2),ncol=n+1)
    for(i in 1:(n+1)){
        for(j in i:(n+1)){
            result1[i,j] <- S0*u^(j-i)*d^(i-1)
        }
    }
    result2[,n+1] <- pmax(0,X-result1[,n+1])
    #倒推求价格
    for(j in n:1){
        for(i in 1:j){
            result2[i,j] <- max(X-result1[i,j],
                                edt*(p*result2[i,j+1]+(1-p)*result2[i+1,j+1]))
        }
    }
    #整理数据
    result <- matrix(rep(" ",(n+1)^2),ncol=n+1)
    for(i in 1:(n+1)){
        for(j in i:(n+1)){
            result[i,j] <- paste0(round(result1[i,j],2),
                                  " || ",round(result2[i,j],2))
        }
    }
    colnames(result) <- 0:n
    return(list(Price=round(result2[1,1],2),
                TABLE=result))
}
##-------------------------------------------------------
#第二张页面布局
page1 <- fluidPage(
    sidebarLayout(
        sidebarPanel(
            titlePanel("输入参数:"),
            # 设置输入参数
            sliderInput('X', '执行价格', 
                        min=0, max=100, value=50, 
                        step=1, round=0),
            sliderInput('S0', '股票现价', 
                        min=0, max=100, value=49, 
                        step=1, round=0),
            sliderInput('r', '无风险年利率(%)', 
                        min=0, max=20, value=5, 
                        step=0.5, round=0),
            sliderInput('sig', '波动率（标准差）', 
                        min=0, max=1, value=0.2, 
                        step=0.05, round=0),
            sliderInput('Ti', '到期时间（周）', 
                        min=10, max=52, value=20, 
                        step=1, round=0),
            sliderInput('dt', '对冲调整间隔（周）', 
                        min=0.1, max=5, value=1, 
                        step=0.1, round=0),
            sliderInput('N', '交易数量（万份）', 
                        min=1, max=50, value=10, 
                        step=1, round=0),
            sliderInput('RT', '模拟重复次数', 
                        min=100, max=2*1e3, value=1000, 
                        step=50, round=0),
        ),
        mainPanel(
            tabsetPanel(
                tabPanel("股价模拟图", plotOutput("plot1")),  
                tabPanel("Delta对冲表", tableOutput("table1")),
                tabPanel("对冲表现分析", tableOutput("re1"))
            )
        )
    )
)

#第3张页面布局
page2 <- fluidPage(
    sidebarLayout(
        sidebarPanel(
            titlePanel("输入参数:"),
            # 设置输入参数
            sliderInput('X1', '执行价格', 
                        min=0, max=100, value=50, 
                        step=1, round=0),
            sliderInput('S01', '股票现价', 
                        min=0, max=100, value=50, 
                        step=1, round=0),
            sliderInput('r1', '无风险年利率(%)', 
                        min=0, max=20, value=10, 
                        step=0.5, round=0),
            sliderInput('sig1', '波动率（标准差）', 
                        min=0, max=1, value=0.4, 
                        step=0.05, round=0),
            sliderInput('Ti1', '到期时间（年）', 
                        min=0.01, max=1, value=0.42, 
                        step=0.01, round=0),
            sliderInput('n1', '二叉树层数', 
                        min=2, max=50, value=5, 
                        step=1, round=0),
            sliderInput('N1', '交易数量（万份）', 
                        min=1, max=50, value=10, 
                        step=1, round=0),
            sliderInput('Maxn1', '模拟最大期数', 
                        min=10, max=100, value=50, 
                        step=1, round=0),
        ),
        mainPanel(
            tabsetPanel(
                tabPanel("二叉树模拟图", tableOutput("Table2")),  
                tabPanel("二叉树定价折线图", plotOutput("plot2")),
                tabPanel("二叉树定价结果分析", tableOutput("re2"))
            )
        )
    )
)
##-------------------------------------------------------

ui <- dashboardPage(
    dashboardHeader(title="《数量金融》大作业"),
    #侧边栏
    dashboardSidebar(
        sidebarMenu(
            menuItem("页面说明", tabName = "desc", icon = icon("dashboard")),
            menuItem("Delta对冲", tabName = "bs", icon = icon("th")),
            menuItem("二叉树定价", tabName = "bitree", icon = icon("file-code-o"))
        )
    ),
    #body
    dashboardBody(
        tabItems(
            tabItem(tabName = "desc",h1("本页面是《数量金融》课程大作业的代码可视化界面。请将本界面与提交的报告一同使用，谢谢！"),
                    h1("编写本页面的全部代码请见：https://github.com/u10010101/Quantitative-Finance/blob/main/app.R"),
                    h1("联系本人：u10010101@163.com"),h1("感谢使用！")),
            tabItem(tabName = "bs",page1),
            tabItem(tabName = "bitree",page2)
        )
    )
    
)

server <- function(input, output){
    
    #页面2
    Delta.result  <-  reactive(DeltaEuroCall(X=as.numeric(input$X),
                                    S0=as.numeric(input$S0),
                                    r=as.numeric(input$r)/100,
                                    sig=as.numeric(input$sig),
                                    Ti=as.numeric(input$Ti),
                                    dt=as.numeric(input$dt),
                                    N=1e4*as.numeric(input$N)))
    output$plot1 <- renderPlot(plot(Delta.result()[[1]],xlab="Time",
                                    ylab="Stock Price",type="o",
                                    xaxt="n",pch=19))
    output$table1 <- renderTable(Delta.result()[[2]])
    output$re1 <- renderTable({
        This.time <- Delta.result()[[3]]
        BS.price <- Black_Scholes_Call(X=as.numeric(input$X),
                                      S0=as.numeric(input$S0),
                                      r=as.numeric(input$r)/100,
                                      sig=as.numeric(input$sig),
                                      Ti=as.numeric(input$Ti),
                                      N=1e4*as.numeric(input$N))
       
        a <- purrr::map_dbl(1:input$RT,~DeltaEuroCall(X=as.numeric(input$X),
                                  S0=as.numeric(input$S0),
                                  r=as.numeric(input$r)/100,
                                  sig=as.numeric(input$sig),
                                  Ti=as.numeric(input$Ti),
                                  dt=as.numeric(input$dt),
                                  N=1e4*as.numeric(input$N))[[3]])
        return(data.frame(
            `此次对冲成本`=This.time,
            `多次平均对冲成本`=mean(a),
            `BS定价价格`=BS.price,
            `多次对冲成本标准差`=sd(a),
            `对冲表现`=sd(a)/BS.price
        ))
    })
    
    #页面3
    bitree.result  <-  reactive(BiAmerCall(X=as.numeric(input$X1),
                                             S0=as.numeric(input$S01),
                                             r=as.numeric(input$r1)/100,
                                             sig=as.numeric(input$sig1),
                                             Ti=as.numeric(input$Ti1),
                                             n=as.numeric(input$n1),
                                             N=as.numeric(input$N1)))
    data3 <- reactive({
        y <- rep(0,input$Maxn1)
        for(k in 1:length(y)){
            y[k] <- BiAmerCall(X=as.numeric(input$X1),
                           S0=as.numeric(input$S01),
                           r=as.numeric(input$r1)/100,
                           sig=as.numeric(input$sig1),
                           Ti=as.numeric(input$Ti1),
                           n=k,
                           N=as.numeric(input$N1))[[1]]
    }
        return(y)})
    output$Table2 <- renderTable({
        re.tmp <- as.data.frame(bitree.result()[[2]])
        ifelse(ncol(re.tmp)>13,
               re.tmp1 <- re.tmp[,c(1:9,(ncol(re.tmp)-3):ncol(re.tmp))],
               re.tmp1 <- re.tmp)
        return(re.tmp1)
    })
    output$plot2 <- renderPlot({
        plot(data3()[-1],type="o",xlab="n",ylab="Option Price")
    })
    output$re2 <- renderTable({
        This.time1 <- bitree.result()[[1]]
        Max.time <- data3()[length(data3())]
        return(data.frame(
            `此次定价结果`=This.time1,
            `最大期数时的定价结果`=Max.time
        ))
    })
    
}


shinyApp(ui, server)
