abstract type Example end

struct Example1 <: Example end
struct Example2 <: Example end

struct Tasks
Ex::Example
end

T1 = Tasks(Example1())
T2 = Tasks(Example2())
function A(x,Task::Example1)
  10*x
end  

function A(x,Task::Example2)
  20*x
end

@show A(5,T1.Ex)
@show A(5,T2.Ex)

