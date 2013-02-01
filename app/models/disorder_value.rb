class DisorderValue
  include DataMapper::Resource
  
  property :disorder_value_id, Serial#, :field => 'disorder_value_id'
  property :disorder_id, Integer, :required => true
  property :aasequence_id, Integer, :required => true
  property :dvalue, Float, :required => true

  #alias :id :disorder_value_id
  #alias :disorder_value_id :id
  belongs_to :disorder, 'Disorder'#, :parent_key => [:disorder_id], :child_key => [:disorder_value_id], :required => true  
  belongs_to :a_asequence, 'AAsequence', :child_key=>[:aasequence_id]
  
  def id
    self.disorder_value_id
  end
end
