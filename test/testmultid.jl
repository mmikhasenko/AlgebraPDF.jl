
#                            _|    _|      _|                    _|  
#  _|_|_|  _|_|    _|    _|  _|  _|_|_|_|                    _|_|_|  
#  _|    _|    _|  _|    _|  _|    _|      _|  _|_|_|_|_|  _|    _|  
#  _|    _|    _|  _|    _|  _|    _|      _|              _|    _|  
#  _|    _|    _|    _|_|_|  _|      _|_|  _|                _|_|_|  


@testset "cross-product PDF" begin
    # test
    pdf1 = pdf((x;p)->x.^2; lims=(-1,2), p=∅)
    pdf2 = pdf((x;p)->x.^4; lims=(-1,2), p=∅)
    X = xProductPDF(;x=pdf1, y=pdf2)
    s = generate(100, X)
    # 
    @test length(s) == 100
    @test hasproperty(s[1], :x)
    @test hasproperty(s[1], :y)
    # 
    s = generate(50, X; Nbins=300)
    @test length(s) == 50
end