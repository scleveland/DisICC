class AAsequence 
  include DataMapper::Resource
  
  property :id, Serial, :field => 'AAsequence_id'
  property :seq_id, Integer, :required => true
  property :amino_acid, String, :length=> 1,  :required => true
  property :original_position, Integer, :required => false
  property :disorder_consensus, Float, :required => false, :default => 0.0
  property :contact_consensus, Float, :required => false, :default => 0.0
  property :contact_positive_consensus, Integer, :required => false, :default => 0.0
  
  alias :AAsequence_id :id
  belongs_to :sequence, 'Sequence', :child_key => [:seq_id]
  #has n, :disorder, 'Disorder', :parent_key=>[:disorder_id]
  has n, :disorder_values, 'DisorderValue', :child_key => [:disorder_value_id]
  has 1, :xdet, 'Xdet', :child_key => :aasequence_id
  has 1, :conseq, 'Conseq', :child_key => :aasequence_id
end
