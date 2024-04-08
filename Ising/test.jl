module Foo
export bar
function bar()
	println("Hello world!")
end
end

using .Foo
bar()
