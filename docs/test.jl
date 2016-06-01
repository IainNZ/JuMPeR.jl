function foo(a::Float64)
    return true
end

function trigger(a, b; c::Float64 = a * b)
    return true
end
trigger(5, 10, c = 3)
#trigger(5, 10)
