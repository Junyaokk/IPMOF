# Online-Retail Inventory Placement and Multi-item Order Fulfillment (IPMOF)

This is our implementation for the paper:
- Online-Retail Inventory Placement and Multi-item Order Fulfillment

We provides the Python codes for data loading, data generation and all algorithms considered in this paper.

If you encounter any significant bugs, then please contact [Junyao Yu](junyaoyu@163.sufe.edu.cn).

## Usage

We have selected two specific retailers' data, namely Retailer 1 and Retailer 2, for evaluating the performance.

To facilitate more accurate estimations of the stochastic program and test different problem sizes for our algorithms, we have developed a scenario generator for the original sets. Define the problem with generated scenarios set as a problem instance, which is served as the input for our algorithms.

The reader can specify the retailer index, replenishment frequency, and the number of scenarios to generate a problem instance.

First, we create an instance base to load and preprocess the data.


```python
from domain.instance import InstanceBase
instance_base = InstanceBase(retailer_id='retailer1', 
                             data_dir='./data/retailer1/', 
                             scenario_unit=7)
```

To optimize the problem instance using the instance base, users need to specify the following parameters to create a task.

- *task_id*: the unique index of the task.
- *instance_base*: prepared data.
- *need_solver*: whether the Gurobi solver should be invoked to compute the final cost based on the solution obtained from our algorithms (True or False).
- *need_opt*: whether the Gurobi solver should be invoked to optimize the problem (True or False).
- *need_benchmark*: whether the product-level benchmark should be invoked to optimize the problem (True or False).
- *capacity_constr*: whether the problem includes a capacity constraint (True or False).
- *capacity_type*: the type of capacity constraint, it can be 'warehouse', 'assortment' or 'mixed'.
- *input_fdc_capacity*: the value of the single capacity constraint.
- *mixed_paras*: the values of the two types capacity constraints.
- *algorithm_paras*: the values of parameters in our algorithms.
- *scenario_num*: the number of scenarios within a problem instance.
- *data_dir*: the directory for storing the task results.
- *repooling_methods_type*: the specified demand aggregation strategy for repooling methods, it can be 'F' or 'T'.

The default parameters for these algorithms are specified in 'default_paras.py'. To introduce the implementation of three algorithms, they are run with the same parameters set. However, this does not guarantee the optimality of the results. For users, we provide a summary of the recommended parameter intervals and values. Further details can be found in Appendix H of our paper.

The relative test functions have been encapsulated, readers can read them in detail. For example, we create three tasks from above instance base for PBS, P\&K and hybrid algorithms.

Here, we simply set the random seed to 3407. Note that under well-tuned parameters setting, the randomness of the problem instances does not affect the performance of our algorithms.


```python
import random
from domain.task import Task
random.seed(3407)
```

- **PH-based Bounded Subgradient (PBS) Algorithm with Repooling Methods**


```python
task_pbs = Task(task_id='TASK_001',
                instance_base=instance_base, 
                capacity_constr=True,
                capacity_type='warehouse',
                need_solver=True,
                need_opt=True,
                need_benchmark=True,
                scenario_num=10,
                input_fdc_capacity=8000)
task_pbs.run(data_dir='./data/retailer1/')
```


```python
task_pbs.run_repooling(repooling_methods_type='F')
task_pbs.run_repooling(repooling_methods_type='T')
```

- **PH-based Algorithm & Knapsack Problem with Dependencies (P&K) with Repooling Methods**


```python
task_pnk = Task(task_id='TASK_002',
                instance_base=instance_base, 
                capacity_constr=True,
                capacity_type='assortment',
                need_solver=True,
                need_opt=True,
                need_benchmark=True,
                scenario_num=10,
                input_fdc_capacity=190
                )
task_pnk.run(data_dir='./data/retailer1/')
```


```python
task_pnk.run_repooling(repooling_methods_type='F')
task_pnk.run_repooling(repooling_methods_type='T')
```

- **Hybrid Algorithm with Repooling Methods**


```python
task_hybrid = Task(task_id='TASK_003',
                instance_base=instance_base, 
                capacity_constr=True,
                capacity_type='mixed',
                need_solver=True,
                need_opt=True,
                need_benchmark=True,
                scenario_num=10,
                mixed_paras={'assortment': 190,
                            'warehouse': 5000}
                )
task_hybrid.run(data_dir='./data/retailer1/')
```


```python
task_hybrid.run_repooling(repooling_methods_type='F')
task_hybrid.run_repooling(repooling_methods_type='T')
```
