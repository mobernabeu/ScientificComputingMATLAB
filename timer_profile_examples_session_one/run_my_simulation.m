function run_my_simulation()

    function do_one_thing(iteration_count)
        % This implementation is a stub that pauses execution for 0.1s.
        pause(0.1)
    end
    function do_another_thing(iteration_count)
        % This implementation is a stub that pauses execution for
        % iteration_count*0.1s.
        pause(0.1*iteration_count)
    end

one_timer = 0;
another_timer = 0;

for iteration_count = 1:10
    t_start=tic;
    do_one_thing(iteration_count);
    one_timer = one_timer + toc(t_start);
    
    t_start=tic;
    do_another_thing(iteration_count);
    another_timer = another_timer + toc(t_start);
end
one_timer
another_timer

end

