class Xdet
  include DataMapper::Resource
  
  property :id, Serial, :field=> 'xdet_id'
  property :aasequence_id, Integer, :required => true
  property :conservation, Float, :required => true
  property :correlation, Float, :required => true

  belongs_to :aasequence, 'AAsequence', :child_key => :aasequence_id
  belongs_to :sequence, :through=>:aasequence
end
