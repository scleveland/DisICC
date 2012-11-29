class DisorderValue
  include DataMapper::Resource
  
  property :id, Serial
  #property :disorder_id, Integer, :required => true
  #property :aasequence_id, Integer, :required => true
  property :dvalue, Float, :required => true

  alias :disorder_value_id :id
  belongs_to :disorder, 'Disorder'#, :parent_key => [:disorder_id], :child_key => [:disorder_value_id], :required => true  
  belongs_to :a_asequence, 'AAsequence'#,  :parent_key => [:AAsequence_id], :child_key=>[:aasequence_id]
end
