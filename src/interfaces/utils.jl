# https://discourse.julialang.org/t/106649
function catch_abbreviated_backtrace()
    bt_caught = catch_backtrace()
    bt_here = backtrace()
    N = length(bt_caught)
    offset = length(bt_here) - N
    while N > 1
        if bt_caught[N] ≢ bt_here[N+offset]
            break
        end
        N = N - 1
    end
    return bt_caught[1:N]
end
