class DisorderValue
  include DataMapper::Resource
  
  property :disorder_value_id, Serial
  property :disorder_id, Integer, :required => true
  property :aasequence_id, Integer, :required => true
  property :dvalue, Float, :required => true

  belongs_to :disorder, 'Disorder', :parent_key => [:disorder_id], :child_key => [:disorder_id]  
  belongs_to :aasequence, 'AAsequence',  :parent_key => [:AAsequence_id], :child_key=>[:aasequence_id]
end
