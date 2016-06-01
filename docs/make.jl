using Documenter, JuMPeR

function Base.Markdown.plain(io::IO, l::Base.Markdown.LaTeX)
    println(io, '$', '$')
    println(io, l.formula)
    println(io, '$', '$')
end

makedocs()
