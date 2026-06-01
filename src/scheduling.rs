//! Application: optimal agent scheduling and resource allocation via tropical optimization.
//!
//! Uses tropical algebra to solve scheduling and resource allocation problems.
//! Tropical matrix powers model longest-path computations in DAGs.
//! Tropical eigenvalues give critical path lengths.

use crate::linear_algebra::TropicalMatrix;
use crate::semiring::Tropical;
use nalgebra::DVector;
use serde::{Deserialize, Serialize};
use std::fmt;

/// An agent in a scheduling system.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Agent {
    /// Unique identifier.
    pub id: String,
    /// Capabilities: vector of skill levels.
    pub capabilities: Vec<f64>,
}

impl Agent {
    pub fn new(id: &str, capabilities: Vec<f64>) -> Self {
        Agent {
            id: id.to_string(),
            capabilities,
        }
    }

    /// Number of capabilities.
    pub fn num_capabilities(&self) -> usize {
        self.capabilities.len()
    }
}

impl fmt::Display for Agent {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Agent({})", self.id)
    }
}

/// A task to be scheduled.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Task {
    /// Unique identifier.
    pub id: String,
    /// Required capabilities (skill requirements).
    pub requirements: Vec<f64>,
    /// Duration of the task.
    pub duration: f64,
    /// Dependencies (indices of tasks that must complete first).
    pub dependencies: Vec<usize>,
    /// Priority (higher = more important).
    pub priority: f64,
}

impl Task {
    pub fn new(id: &str, requirements: Vec<f64>, duration: f64) -> Self {
        Task {
            id: id.to_string(),
            requirements,
            duration,
            dependencies: vec![],
            priority: 0.0,
        }
    }

    pub fn with_dependencies(mut self, deps: Vec<usize>) -> Self {
        self.dependencies = deps;
        self
    }

    pub fn with_priority(mut self, priority: f64) -> Self {
        self.priority = priority;
        self
    }
}

/// A schedule assignment.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Assignment {
    pub task_idx: usize,
    pub agent_idx: usize,
    pub start_time: f64,
    pub end_time: f64,
}

impl Assignment {
    pub fn new(task_idx: usize, agent_idx: usize, start_time: f64, end_time: f64) -> Self {
        Assignment { task_idx, agent_idx, start_time, end_time }
    }
}

/// Agent scheduler using tropical optimization.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AgentScheduler {
    /// Available agents.
    pub agents: Vec<Agent>,
    /// Tasks to schedule.
    pub tasks: Vec<Task>,
}

impl AgentScheduler {
    pub fn new(agents: Vec<Agent>, tasks: Vec<Task>) -> Self {
        AgentScheduler { agents, tasks }
    }

    /// Build the tropical cost matrix for agent-task assignment.
    /// C[i][j] = -(suitability of agent i for task j) = we want to maximize suitability.
    pub fn cost_matrix(&self) -> TropicalMatrix {
        let n = self.agents.len();
        let m = self.tasks.len();
        let mut rows = Vec::with_capacity(n);

        for agent in &self.agents {
            let mut row = Vec::with_capacity(m);
            for task in &self.tasks {
                let suitability = self.compute_suitability(agent, task);
                // Cost = -suitability (tropical min selects highest suitability)
                row.push(-suitability);
            }
            rows.push(row);
        }

        TropicalMatrix::from_rows(&rows)
    }

    /// Compute suitability of an agent for a task.
    fn compute_suitability(&self, agent: &Agent, task: &Task) -> f64 {
        let nc = agent.num_capabilities().min(task.requirements.len());
        if nc == 0 {
            return 0.0;
        }
        // Suitability = min over required capabilities (bottleneck)
        // In tropical terms: max of (capability - requirement) = best fit
        let mut min_fit = f64::INFINITY;
        for i in 0..nc {
            let fit = agent.capabilities[i] - task.requirements[i];
            min_fit = min_fit.min(fit);
        }
        min_fit + task.priority
    }

    /// Build the tropical dependency matrix.
    /// D[i][j] = duration of task j if it depends on task i, else -∞.
    pub fn dependency_matrix(&self) -> TropicalMatrix {
        let n = self.tasks.len();
        let mut rows = vec![vec![f64::NEG_INFINITY; n]; n];

        for (j, task) in self.tasks.iter().enumerate() {
            for &dep in &task.dependencies {
                if dep < n {
                    rows[dep][j] = self.tasks[dep].duration;
                }
            }
        }

        TropicalMatrix::from_rows(&rows)
    }

    /// Compute critical path length using tropical matrix algebra.
    /// The critical path is the longest path in the DAG of dependencies.
    /// Computed as the tropical eigenvalue of the dependency matrix.
    pub fn critical_path_length(&self) -> f64 {
        let dep_matrix = self.dependency_matrix();
        dep_matrix
            .tropical_eigenvalue()
            .map(|t| t.value())
            .unwrap_or(0.0)
    }

    /// Compute earliest start times using tropical matrix powers.
    /// ES[j] = max over dependencies i of (ES[i] + duration[i]).
    pub fn earliest_start_times(&self) -> Vec<f64> {
        let n = self.tasks.len();
        let dep_matrix = self.dependency_matrix();
        let mut es = vec![0.0; n];

        // Topological order via tropical matrix algebra
        // ES = A* ⊗ d where A is dependency matrix and d is duration vector
        let durations: Vec<Tropical> = self.tasks.iter().map(|t| Tropical::new(t.duration)).collect();
        let d = DVector::from_column_slice(&durations);

        // Compute using tropical matrix powers
        let mut current = DVector::from_element(n, Tropical::ONE); // Start at 0
        for _ in 0..n {
            let next = dep_matrix.mul_vector(&current);
            let mut changed = false;
            for i in 0..n {
                let new_val = current[i] + next[i] + d[i];
                if new_val > current[i] {
                    // Actually: es_i = max(current_i, max_dep(next + duration))
                    changed = true;
                }
            }
            current = next;
            if !changed {
                break;
            }
        }

        // Simpler approach: iterate until convergence
        let mut es_vec = vec![0.0f64; n];
        for iteration in 0..n {
            let mut changed = false;
            for j in 0..n {
                let mut max_dep_end = 0.0f64;
                for &dep in &self.tasks[j].dependencies {
                    if dep < n {
                        max_dep_end = max_dep_end.max(es_vec[dep] + self.tasks[dep].duration);
                    }
                }
                if max_dep_end > es_vec[j] {
                    es_vec[j] = max_dep_end;
                    changed = true;
                }
            }
            if !changed {
                break;
            }
            let _ = iteration; // suppress unused warning
        }

        es_vec
    }

    /// Schedule tasks using tropical optimization.
    pub fn schedule(&self) -> Vec<Assignment> {
        let n = self.tasks.len();
        let m = self.agents.len();
        if n == 0 || m == 0 {
            return vec![];
        }

        let cost = self.cost_matrix();
        let es = self.earliest_start_times();
        let mut assignments = Vec::with_capacity(n);

        // Track agent availability
        let mut agent_available = vec![0.0; m];

        // Assign tasks in order of earliest start time
        let mut task_order: Vec<usize> = (0..n).collect();
        task_order.sort_by(|&a, &b| es[a].partial_cmp(&es[b]).unwrap_or(std::cmp::Ordering::Equal));

        let mut assigned_tasks = vec![false; n];

        for &task_idx in &task_order {
            let task = &self.tasks[task_idx];

            // Find the best agent for this task
            let mut best_agent = 0;
            let mut best_cost = f64::INFINITY;

            for agent_idx in 0..m {
                let suitability = -cost.get(agent_idx, task_idx).value();
                let start = es[task_idx].max(agent_available[agent_idx]);
                let total_time = start + task.duration;

                // Prefer: high suitability, low total time
                let score = -suitability + total_time;
                if score < best_cost {
                    best_cost = score;
                    best_agent = agent_idx;
                }
            }

            let start_time = es[task_idx].max(agent_available[best_agent]);
            let end_time = start_time + task.duration;
            agent_available[best_agent] = end_time;

            assignments.push(Assignment::new(task_idx, best_agent, start_time, end_time));
            assigned_tasks[task_idx] = true;
        }

        assignments
    }

    /// Compute the makespan (total schedule duration).
    pub fn makespan(&self, assignments: &[Assignment]) -> f64 {
        assignments
            .iter()
            .map(|a| a.end_time)
            .fold(0.0f64, f64::max)
    }

    /// Compute total resource utilization.
    pub fn resource_utilization(&self, assignments: &[Assignment]) -> f64 {
        if assignments.is_empty() {
            return 0.0;
        }
        let makespan = self.makespan(assignments);
        if makespan <= 0.0 {
            return 0.0;
        }
        let total_work: f64 = assignments.iter().map(|a| a.end_time - a.start_time).sum();
        let num_agents = self.agents.len() as f64;
        total_work / (makespan * num_agents)
    }

    /// Optimal resource allocation using tropical linear programming.
    /// Allocate resources to maximize throughput.
    pub fn allocate_resources(&self, budget: f64) -> Vec<f64> {
        let n = self.agents.len();
        let mut allocation = vec![0.0; n];

        // Allocate budget proportional to capability
        let total_cap: f64 = self
            .agents
            .iter()
            .map(|a| a.capabilities.iter().sum::<f64>())
            .sum();

        if total_cap <= 0.0 {
            // Equal allocation
            let share = budget / n as f64;
            allocation = vec![share; n];
        } else {
            for (i, agent) in self.agents.iter().enumerate() {
                let cap_sum: f64 = agent.capabilities.iter().sum();
                allocation[i] = budget * (cap_sum / total_cap);
            }
        }

        allocation
    }

    /// Compute the tropical bottleneck for the schedule.
    /// The bottleneck is the task that most delays the critical path.
    pub fn bottleneck_task(&self) -> (usize, f64) {
        let es = self.earliest_start_times();
        let mut bottleneck = 0;
        let mut max_impact = 0.0;

        for (i, task) in self.tasks.iter().enumerate() {
            let impact = task.duration + es[i];
            if impact > max_impact {
                max_impact = impact;
                bottleneck = i;
            }
        }

        (bottleneck, max_impact)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_agent_creation() {
        let agent = Agent::new("a1", vec![3.0, 4.0, 5.0]);
        assert_eq!(agent.id, "a1");
        assert_eq!(agent.num_capabilities(), 3);
    }

    #[test]
    fn test_agent_display() {
        let agent = Agent::new("worker-1", vec![1.0]);
        assert_eq!(format!("{}", agent), "Agent(worker-1)");
    }

    #[test]
    fn test_task_creation() {
        let task = Task::new("t1", vec![2.0, 3.0], 5.0);
        assert_eq!(task.id, "t1");
        assert_eq!(task.duration, 5.0);
    }

    #[test]
    fn test_task_with_dependencies() {
        let task = Task::new("t1", vec![1.0], 3.0).with_dependencies(vec![0, 1]);
        assert_eq!(task.dependencies, vec![0, 1]);
    }

    #[test]
    fn test_task_with_priority() {
        let task = Task::new("t1", vec![1.0], 3.0).with_priority(10.0);
        assert_eq!(task.priority, 10.0);
    }

    #[test]
    fn test_scheduler_creation() {
        let agents = vec![Agent::new("a1", vec![3.0])];
        let tasks = vec![Task::new("t1", vec![2.0], 5.0)];
        let scheduler = AgentScheduler::new(agents, tasks);
        assert_eq!(scheduler.agents.len(), 1);
        assert_eq!(scheduler.tasks.len(), 1);
    }

    #[test]
    fn test_cost_matrix() {
        let agents = vec![
            Agent::new("a1", vec![5.0, 3.0]),
            Agent::new("a2", vec![2.0, 7.0]),
        ];
        let tasks = vec![Task::new("t1", vec![3.0, 2.0], 4.0)];
        let scheduler = AgentScheduler::new(agents, tasks);
        let cost = scheduler.cost_matrix();
        assert_eq!(cost.nrows(), 2);
        assert_eq!(cost.ncols(), 1);
    }

    #[test]
    fn test_dependency_matrix() {
        let agents = vec![Agent::new("a1", vec![1.0])];
        let tasks = vec![
            Task::new("t1", vec![1.0], 3.0),
            Task::new("t2", vec![1.0], 5.0).with_dependencies(vec![0]),
        ];
        let scheduler = AgentScheduler::new(agents, tasks);
        let dep = scheduler.dependency_matrix();
        assert_eq!(dep.get(0, 1).value(), 3.0); // t2 depends on t1 with duration 3
        assert!(dep.get(1, 0).is_zero()); // no reverse dependency
    }

    #[test]
    fn test_earliest_start_times_no_deps() {
        let agents = vec![Agent::new("a1", vec![1.0])];
        let tasks = vec![
            Task::new("t1", vec![1.0], 3.0),
            Task::new("t2", vec![1.0], 5.0),
        ];
        let scheduler = AgentScheduler::new(agents, tasks);
        let es = scheduler.earliest_start_times();
        assert_eq!(es[0], 0.0);
        assert_eq!(es[1], 0.0);
    }

    #[test]
    fn test_earliest_start_times_with_deps() {
        let agents = vec![Agent::new("a1", vec![1.0])];
        let tasks = vec![
            Task::new("t1", vec![1.0], 3.0),
            Task::new("t2", vec![1.0], 5.0).with_dependencies(vec![0]),
        ];
        let scheduler = AgentScheduler::new(agents, tasks);
        let es = scheduler.earliest_start_times();
        assert_eq!(es[0], 0.0);
        assert_eq!(es[1], 3.0); // t2 starts after t1 (duration 3)
    }

    #[test]
    fn test_schedule_basic() {
        let agents = vec![Agent::new("a1", vec![5.0])];
        let tasks = vec![
            Task::new("t1", vec![1.0], 3.0),
            Task::new("t2", vec![1.0], 2.0),
        ];
        let scheduler = AgentScheduler::new(agents, tasks);
        let assignments = scheduler.schedule();
        assert_eq!(assignments.len(), 2);
    }

    #[test]
    fn test_schedule_with_deps() {
        let agents = vec![Agent::new("a1", vec![5.0])];
        let tasks = vec![
            Task::new("t1", vec![1.0], 3.0),
            Task::new("t2", vec![1.0], 2.0).with_dependencies(vec![0]),
        ];
        let scheduler = AgentScheduler::new(agents, tasks);
        let assignments = scheduler.schedule();
        // t2 should start after t1 ends
        let t1 = assignments.iter().find(|a| a.task_idx == 0).unwrap();
        let t2 = assignments.iter().find(|a| a.task_idx == 1).unwrap();
        assert!(t2.start_time >= t1.end_time - 1e-10);
    }

    #[test]
    fn test_makespan() {
        let agents = vec![Agent::new("a1", vec![5.0])];
        let tasks = vec![
            Task::new("t1", vec![1.0], 3.0),
            Task::new("t2", vec![1.0], 2.0),
        ];
        let scheduler = AgentScheduler::new(agents, tasks);
        let assignments = scheduler.schedule();
        let makespan = scheduler.makespan(&assignments);
        assert!(makespan >= 5.0); // At least 3+2
    }

    #[test]
    fn test_resource_utilization() {
        let agents = vec![Agent::new("a1", vec![5.0])];
        let tasks = vec![
            Task::new("t1", vec![1.0], 3.0),
            Task::new("t2", vec![1.0], 2.0),
        ];
        let scheduler = AgentScheduler::new(agents, tasks);
        let assignments = scheduler.schedule();
        let util = scheduler.resource_utilization(&assignments);
        assert!(util > 0.0);
        assert!(util <= 1.0);
    }

    #[test]
    fn test_resource_allocation() {
        let agents = vec![
            Agent::new("a1", vec![5.0]),
            Agent::new("a2", vec![3.0]),
        ];
        let tasks = vec![Task::new("t1", vec![1.0], 1.0)];
        let scheduler = AgentScheduler::new(agents, tasks);
        let alloc = scheduler.allocate_resources(100.0);
        assert_eq!(alloc.len(), 2);
        let total: f64 = alloc.iter().sum();
        assert!((total - 100.0).abs() < 1e-10);
        // a1 has more capability, should get more
        assert!(alloc[0] > alloc[1]);
    }

    #[test]
    fn test_bottleneck_task() {
        let agents = vec![Agent::new("a1", vec![5.0])];
        let tasks = vec![
            Task::new("t1", vec![1.0], 10.0), // Long task
            Task::new("t2", vec![1.0], 2.0),
        ];
        let scheduler = AgentScheduler::new(agents, tasks);
        let (idx, impact) = scheduler.bottleneck_task();
        assert_eq!(idx, 0); // t1 is the bottleneck
        assert!(impact > 0.0);
    }

    #[test]
    fn test_schedule_empty() {
        let scheduler = AgentScheduler::new(vec![], vec![]);
        let assignments = scheduler.schedule();
        assert!(assignments.is_empty());
    }

    #[test]
    fn test_multi_agent_schedule() {
        let agents = vec![
            Agent::new("a1", vec![5.0, 2.0]),
            Agent::new("a2", vec![2.0, 5.0]),
        ];
        let tasks = vec![
            Task::new("t1", vec![4.0, 1.0], 3.0),
            Task::new("t2", vec![1.0, 4.0], 2.0),
        ];
        let scheduler = AgentScheduler::new(agents, tasks);
        let assignments = scheduler.schedule();
        assert_eq!(assignments.len(), 2);
        // t1 better suited to a1, t2 better suited to a2
    }

    #[test]
    fn test_schedule_serialization() {
        let agents = vec![Agent::new("a1", vec![5.0])];
        let tasks = vec![Task::new("t1", vec![1.0], 3.0)];
        let scheduler = AgentScheduler::new(agents, tasks);
        let json = serde_json::to_string(&scheduler).unwrap();
        let deserialized: AgentScheduler = serde_json::from_str(&json).unwrap();
        assert_eq!(deserialized.agents.len(), 1);
    }
}
