# q_learning

How do humans act in uncertain environments? What are the computational principles that might guide such decision making? 
And in more neurobiological terms, what are the neural structures involved in computing and assessing reward and risk?
For the purposes of this study I constructed a two-dimensional maze for the agent to navigate in (see figure). 
The starting position is at coordinate (5, 6) marked in yellow. The purpose of the agent is to exit 
(6, 1) the maze via a route that yields an optimal payoff. 
By optimal here I mean the minimum punishment/maximum reward possible.

The agent has a choice of routes, route A and route B. the former leads to the exit via a long detour whilst the latter leads to 
the exit via a much shorter (but riskier) route. 
However there is a trap along route B which forces the agent to make decision based on its prior knowledge of the incidence 
of the traps. The occurrence of the trap within the maze was modelled through a Poisson process. By varying 
the parameters of the Poisson process I was able to set the risk involved in taking this specific route which led to the exit. 

There were thus two independent variables: the risk setting of the trap (as quantified by the λ value) 
and the degree of risk seeking undertaken by the subject (as quantified by either ε for ε-greedy or α for the SOFTMAX algorithm). 
The learning constant, α, was set at 1 for all {λ, ε|α} pairs.

The following observable variables were collected for analysis:

1.	Total payoff
2.	Duration per episode
3.	Prediction error, that is: |Q(s,a)t – Q(s,a)t-1|

These variables were then used in order to gain a better understanding of how the agent is forming 
an accurate model of the statistics of its world.
