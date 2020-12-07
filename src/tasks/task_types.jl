# Generate Task and TaskGDS types with names e.g. Attractor and AttractorGDS
# and no other fields
macro standard_task_type(name)
    esc(quote
        struct $name{F} <: Task{F}
            task_map::F
        end

        struct $(Symbol("$name",:GDS)){F} <: TaskGDS{F}
            task_map::F
        end
    end)
end

@standard_task_type BlankTask
@standard_task_type Damping
@standard_task_type Attractor
@standard_task_type CollisionAvoidance
@standard_task_type VelocityLimit

macro ext_standard_task_type(name)
    esc(quote
        struct $name{F} <: PBDS.Task{F}
            task_map::F
        end

        struct $(Symbol("$name",:GDS)){F} <: PBDS.TaskGDS{F}
            task_map::F
        end
    end)
end