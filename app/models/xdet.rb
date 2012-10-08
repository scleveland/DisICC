class Xdet
  include DataMapper::Resource
  
  property :xdet_id, Serial
  property :aasequence_id, Integer, :required => true
  property :conservation, Float, :required => true
  property :correlation, Float, :required => true

  

end
