# FtsZ_acc
The repository contains simulation codes used in the paper "Determining the rate-limiting processes for cell division in Escherichia coli", Nature Communications (2024).

Three executable MATLAB files and Get_params folder for 
concurrent_prod_rate_change_timer_from_init.m - Concurrent processes model is simulated over a population of cells. The concurrent process is at the constriction. Constriction is controlled by a timer from initiation and accumulation of an division protein from birth. Division happens after a constant time from constriction. Growth rate is different before and after constriction. At certain time, there is over expression of the division protein. The production rate of division protein changes and the constriction time changes but the the threshold to which the division protein is accumulated at constriction remains constant.

concurrent_prod_rate_change_2_div_protein.m - Concurrent processes model is simulated over a population of cells. The concurrent process is at the constriction. Constriction is controlled by a two division proteins being accumulated from birth. Division happens after a constant time from constriction. Growth rate is different before and after constriction. At certain time, there is over expression of one of the division proteins. The production rate of division protein changes and the constriction time changes but the threshold to which the division protein is accumulated at constriction remains constant.

concurrent_prod_rate_change_2_div_protein.m - Concurrent processes model is simulated over a population of cells. The concurrent process is at the constriction. Constriction is controlled by a timer from initiation and accumulation of an division protein from birth. Division happens after a constant time from constriction. Growth rate is different before and after constriction. At certain time, the threshold to which the division protein accumulates changes.

Get_params - Run the call.m function to get the parameters of two simulation models. Finds the best parameters for the threshold level at constriction of FtsZ and time elapsed between initiation and constriction for the concurrent_prod_rate_change_timer_from_init.m program. Also finds the best parameter for the threshold level of the other division protein at constriction in the concurrent_prod_rate_change_2_div_protein.m program.

