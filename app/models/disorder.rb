class Disorder
  include DataMapper::Resource
  
  property :id, Serial, :field=> 'disorder_id'
  #property :seq_id, Integer, :required => true
  property :disorder_type, String, :required => true
  property :version, Integer, :required => false
  property :seq_id, Integer
  property :deleted_at, ParanoidDateTime
  
  alias :disorder_id :id
  belongs_to :sequence, 'Sequence', :child_key =>[:seq_id], :parent_key => [:seq_id]
  has n, :disorder_values, 'DisorderValue'
end
