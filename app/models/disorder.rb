class Disorder
  include DataMapper::Resource
  
  property :id, Serial
  #property :seq_id, Integer, :required => true
  property :disorder_type, String, :required => true
  property :version, Integer, :required => false
  property :seq_id, Integer
  
  alias :disorder_id :id
  belongs_to :sequence, 'Sequence'
  has n, :disorder_values, 'DisorderValue'
end
